"""Shared helpers for the prepare-ncbi-submission commands."""

import logging
import re
import shutil
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

# NCBI limits submissions to 3 000 sequences per file.
NCBI_BATCH_SIZE = 2999

PLAIN_HEADER_VIRUSES = frozenset(
    {
        "SARS-CoV-2",
        "Dengue virus type 1",
        "Dengue virus type 2",
        "Dengue virus type 3",
        "Dengue virus type 4",
        "Influenza A H1N1",
        "Influenza A H3N2",
        "Influenza B virus (B/Lee/1940)",
        "Influenza C virus (C/Ann Arbor/1/50)",
        "Norovirus GI",
        "Norovirus GII",
        "Norovirus GIII",
        "Norovirus GIV",
        "Norovirus GV",
        "Norovirus GVI",
    }
)


_INFLUENZA_RE = re.compile(r"^Influenza\s+(A|B|C)", re.IGNORECASE)
_DENGUE_RE = re.compile(r"Dengue virus type\s+(\d)", re.IGNORECASE)
_NOROVIRUS_RE = re.compile(
    r"^Norovirus\b.*?(G(?:VI|V|IV|III|II|I))\b", re.IGNORECASE
)
# Characters forbidden in FASTA headers / filenames
_UNSAFE_RE = re.compile(r"[^\x00-\x7F]|[|/\\\s]")


def sanitize_seq_name(name: str) -> tuple[str, bool]:
    """
    Sanitize a sequence name by replacing unsafe characters.

    Replaces non-ASCII chars, slashes, pipes, and whitespace with
    underscores. Preserves trailing [Organism=...] tags.
    Returns both the cleaned name and a flag indicating
    whether any substitutions were made.

    Args:
        name: Original sequence name.

    Returns:
        Tuple of (sanitized_name, was_changed).
    """
    m = re.search(r"(\s+\[Organism=.*?\])$", name)
    suffix = ""
    base_name = name
    if m:
        suffix = m.group(1)
        base_name = name[: m.start()]

    clean = _UNSAFE_RE.sub("_", base_name)
    final_clean = clean + suffix
    return final_clean, final_clean != name


def load_results(results_path: Path) -> pd.DataFrame:
    """
    Load a ViralQC results file (TSV, CSV or JSON).

    Args:
        results_path: Path to the results file.

    Returns:
        DataFrame with all result columns.

    Raises:
        ValueError: If the file extension is not recognised.
    """
    suffix = results_path.suffix.lower()
    if suffix == ".tsv":
        return pd.read_csv(results_path, sep="\t", dtype=str, low_memory=False)
    elif suffix == ".csv":
        return pd.read_csv(results_path, sep=";", dtype=str, low_memory=False)
    elif suffix == ".json":
        df = pd.read_json(results_path, orient="table")
        return df.astype(str)
    else:
        raise ValueError(
            f"Unsupported results format '{suffix}'. Use .tsv, .csv, or .json."
        )


def load_fasta(fasta_path: Path) -> dict:
    """
    Parse a FASTA file into a dictionary keyed by sequence ID.

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        Dictionary mapping sequence ID (str) to SeqRecord.
    """
    return SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))


def _n_fraction(seq: str) -> float:
    """Return the fraction of N/n characters in a sequence string."""
    upper = seq.upper()
    return upper.count("N") / len(upper) if upper else 1.0


def filter_sequences(
    fasta_records: dict,
    results_df: pd.DataFrame,
    min_length: int = 150,
    max_n_fraction: float = 0.5,
) -> tuple[list, list]:
    """
    Filter sequences by length and N-content.

    Args:
        fasta_records: Dict of seqName → SeqRecord from the target FASTA.
        results_df: ViralQC results DataFrame (used to cross-reference names).
        min_length: Minimum acceptable sequence length (inclusive).
        max_n_fraction: Maximum acceptable fraction of N bases (exclusive).

    Returns:
        Tuple of (kept_records, dropped_records) where each element is a list
        of dicts with keys 'seqName', 'reason', and optionally 'length'/'n_pct'.
    """
    kept: list[SeqRecord] = []
    dropped: list[dict] = []

    for seq_id, record in fasta_records.items():
        seq_str = str(record.seq)
        length = len(seq_str)
        n_frac = _n_fraction(seq_str)

        if length < min_length:
            dropped.append(
                {
                    "seqName": seq_id,
                    "length": length,
                    "n_pct": round(n_frac * 100, 2),
                    "reason": f"Sequence shorter than {min_length} nt ({length} nt)",
                }
            )
        elif n_frac >= max_n_fraction:
            dropped.append(
                {
                    "seqName": seq_id,
                    "length": length,
                    "n_pct": round(n_frac * 100, 2),
                    "reason": f"N content ≥ {int(max_n_fraction * 100)}% ({n_frac*100:.1f}%)",
                }
            )
        else:
            kept.append(record)

    return kept, dropped


def _batch_path(base_path: Path, batch_index: int) -> Path:
    """
    Return a numbered batch path, e.g. ``sequences.fasta → sequences.1.fasta``.

    Args:
        base_path: The original output path (e.g. ``out_dir / "sequences.fasta"``).
        batch_index: 1-based batch number.

    Returns:
        Path with the batch number inserted before the extension.
    """
    return base_path.with_suffix(f".{batch_index}{base_path.suffix}")


def _write_fasta_single(
    records: list,
    out_path: Path,
    header_fn=None,
    rename_log: list | None = None,
) -> None:
    """Write a single batch of SeqRecords to a FASTA file."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="ascii", errors="replace") as fh:
        for rec in records:
            raw_header = header_fn(rec) if header_fn else rec.id
            clean_header, changed = sanitize_seq_name(raw_header)
            if changed and rename_log is not None:
                rename_log.append(
                    {
                        "original": raw_header,
                        "sanitized": clean_header,
                    }
                )
            fh.write(f">{clean_header}\n{rec.seq}\n")


def write_fasta(
    records: list,
    out_path: Path,
    header_fn=None,
    rename_log: list | None = None,
) -> None:
    """
    Write SeqRecords to FASTA file(s), optionally applying a custom header.

    When the number of records exceeds :data:`NCBI_BATCH_SIZE` (3 000),
    the output is automatically split into numbered files
    (e.g. ``sequences.1.fasta``, ``sequences.2.fasta``, …).

    Non-ASCII / unsafe characters in the final header are sanitized and
    logged to *rename_log* (a list of dicts) if provided.

    Args:
        records: List of SeqRecord objects.
        out_path: Destination file path.
        header_fn: Optional callable (record) → str to compute the header.
            When None, the record.id is used.
        rename_log: Optional list to append rename events to.
    """
    if len(records) <= NCBI_BATCH_SIZE:
        _write_fasta_single(records, out_path, header_fn, rename_log)
        return

    for i in range(0, len(records), NCBI_BATCH_SIZE):
        batch = records[i : i + NCBI_BATCH_SIZE]
        batch_idx = i // NCBI_BATCH_SIZE + 1
        _write_fasta_single(batch, _batch_path(out_path, batch_idx), header_fn, rename_log)
        logger.info(
            "Batch %d: %d sequences written to %s",
            batch_idx,
            len(batch),
            _batch_path(out_path, batch_idx),
        )


def write_dropped_table(dropped: list, out_path: Path) -> None:
    """
    Write the table of dropped sequences to a TSV file.

    Args:
        dropped: List of dicts with at least 'seqName' and 'reason'.
        out_path: Destination file path.
    """
    if not dropped:
        if out_path.exists():
            out_path.unlink()
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(dropped).to_csv(out_path, sep="\t", index=False)


def copy_tbl(tbl_path_str: str, dest_dir: Path, dest_name: str | None = None) -> bool:
    """
    Copy a TBL file into a destination directory.

    Args:
        tbl_path_str: String path to the TBL file (may be empty/NaN).
        dest_dir: Directory to copy the file into.
        dest_name: Optional new filename for the copied TBL file.

    Returns:
        True if the file was copied successfully, False otherwise.
    """
    if not tbl_path_str or pd.isna(tbl_path_str) or str(tbl_path_str).strip() == "":
        return False
    src = Path(str(tbl_path_str).strip())
    if not src.exists():
        logger.warning("TBL file not found: %s", src)
        return False
    dest_dir.mkdir(parents=True, exist_ok=True)
    target_name = dest_name if dest_name else src.name
    shutil.copy2(src, dest_dir / target_name)
    return True


def write_rename_log(rename_log: list, out_path: Path) -> None:
    """
    Write header rename events to a TSV log file.

    Args:
        rename_log: List of dicts with 'original' and 'sanitized' keys.
        out_path: Destination log file path.
    """
    if not rename_log:
        if out_path.exists():
            out_path.unlink()
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rename_log).to_csv(out_path, sep="\t", index=False)


def is_plain_header_virus(virus_name: str) -> bool:
    """
    Return True if the virus uses a plain seqName header (no [Organism=…]).

    Covers SARS-CoV-2, all Dengue types, and all Influenza A/B/C.

    Args:
        virus_name: Value from the 'virus' column of the results file.

    Returns:
        True if the virus should not include an Organism qualifier.
    """
    if virus_name in PLAIN_HEADER_VIRUSES:
        return True
    if _INFLUENZA_RE.match(virus_name):
        return True
    if _NOROVIRUS_RE.match(virus_name):
        return True
    return False


_METADATA_REQUIRED = {
    "Sequence_ID",
    "geo_loc_name",
    "host",
    "isolate",
    "collection-date",
    "isolation-source",
}

_STANDARD_OUTPUT_COLS = [
    "Sequence_ID",
    "geo_loc_name",
    "host",
    "isolate",
    "collection-date",
    "isolation-source",
    "serotype",
]
# SARS-CoV-2 submissions don't include isolation-source or serotype
_SC2_OUTPUT_COLS = [
    "Sequence_ID",
    "geo_loc_name",
    "host",
    "isolate",
    "collection-date",
]
_NOROVIRUS_OUTPUT_COLS = [
    "Sequence_ID",
    "collection-date",
    "genotype",
    "geo_loc_name",
    "host",
    "isolate",
    "isolation-source",
]
_CUSTOM_RENAME_MAP = {
    "geo_loc_name": "Country (geo_loc_name)",
    "host": "Host",
    "isolate": "Isolate",
    "collection-date": "Collection_date",
    "isolation-source": "Isolation_source",
}
_CUSTOM_OUTPUT_COLS = [
    "Sequence_ID",
    "Country (geo_loc_name)",
    "Host",
    "Isolate",
    "Collection_date",
    "Isolation_source",
]


def load_metadata(path: "Path", is_standard: bool = True) -> "pd.DataFrame":
    """
    Load and validate a sample metadata CSV file.

    For standard viruses, enforces all required columns.
    For custom viruses, only enforces `Sequence_ID`.
    Also enforces that `Sequence_ID` values are under 25 characters.

    Args:
        path: Path to the CSV file.
        is_standard: True if checking for standard virus requirements.

    Returns:
        Validated DataFrame (all columns as str).

    Raises:
        ValueError: If any required column is missing, or if any Sequence_ID is >= 25 chars.
    """
    df = pd.read_csv(path, dtype=str)
    required = _METADATA_REQUIRED if is_standard else {"Sequence_ID"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Metadata file is missing required columns: {sorted(missing)}\n"
            f"Found: {sorted(df.columns)}"
        )

    # Check Sequence_ID length, ncbi rule
    long_ids = df[df["Sequence_ID"].str.len() >= 25]["Sequence_ID"].tolist()
    if long_ids:
        raise ValueError(
            f"Sequence_ID values must be less than 25 characters. "
            f"Found {len(long_ids)} invalid IDs, e.g.: {long_ids[:5]}"
        )

    return df


def _serotype_fields(virus: str, clade: str) -> dict:
    """
    Derive serotype (and optional genotype) fields for a standard virus row.

    - SARS-CoV-2  → ``{serotype: clade}``
    - Dengue       → ``{serotype: "1"/"2"/"3"/"4", genotype: clade}``
    - Influenza A  → ``{serotype: "H<n>N<n>"}``; empty for B/C.

    Args:
        virus: Value from the ``virus`` column.
        clade: Value from the ``clade`` column.

    Returns:
        Dict of extra columns to merge into the output row.
    """
    clade_clean = str(clade) if clade and str(clade) not in ("nan", "") else ""
    if "SARS-CoV-2" in virus:
        return {}
    m_dengue = _DENGUE_RE.search(virus)
    if m_dengue:
        return {"genotype": m_dengue.group(1), "serotype": clade_clean}
    m_noro = _NOROVIRUS_RE.search(virus)
    if m_noro:
        # Extract the genotype (e.g. GII.17, GVI.1, GI) from the virus string
        noro_after = re.sub(r"^Norovirus\s+", "", virus, flags=re.IGNORECASE)
        geno_match = re.search(
            r"(G(?:VI|V|IV|III|II|I)(?:\.\w+)?)", noro_after, re.IGNORECASE
        )
        noro_genotype = geno_match.group(1) if geno_match else m_noro.group(1)
        return {"genotype": noro_genotype}
    m_flu = re.search(r"(H\d+N\d+)", virus)
    if m_flu:
        return {"serotype": m_flu.group(1)}
    return {"serotype": ""}


def write_submission_metadata(
    kept_records: list,
    results_df: "pd.DataFrame",
    metadata_df: "pd.DataFrame | None",
    out_path: "Path",
    is_standard: bool = True,
    organism: str = "",
    base_cols: "list | None" = None,
) -> None:
    """
    Build and write the NCBI submission metadata CSV alongside sequences.fasta.

    For *standard* viruses (SARS-CoV-2, Dengue, Influenza) the output columns
    are: ``Sequence_ID, geo_loc_name, host, isolate, collection-date,
    isolation-source, serotype``.

    For *custom* viruses the output columns are: ``Sequence_ID,
    Country (geo_loc_name), Host, Isolate, Collection_date, Isolation_source``.

    Args:
        kept_records: SeqRecord list of sequences that passed filtering.
        results_df: Full ViralQC results DataFrame.
        metadata_df: User-supplied metadata DataFrame (may be None if optional).
        out_path: Destination CSV path.
        is_standard: True for standard-virus format, False for custom.
        organism: Organism string to fill the Organism column (custom only).
        base_cols: Optional column list to override the default.
    """
    if not kept_records:
        return

    seq_ids = [r.id for r in kept_records]

    if metadata_df is None:
        meta = pd.DataFrame({"Sequence_ID": seq_ids})
    else:
        meta = metadata_df[metadata_df["Sequence_ID"].isin(seq_ids)].copy()

        if meta.empty:
            logger.warning(
                "No metadata found for the kept sequences — metadata CSV not written."
            )
            return

    if is_standard:
        results_sub = results_df[results_df["seqName"].isin(seq_ids)][
            ["seqName", "virus", "clade"]
        ].copy()
        results_sub = results_sub.rename(columns={"seqName": "Sequence_ID"})
        meta = meta.merge(results_sub, on="Sequence_ID", how="left")

        extra = meta.apply(
            lambda r: _serotype_fields(
                str(r.get("virus", "") or ""),
                str(r.get("clade", "") or ""),
            ),
            axis=1,
            result_type="expand",
        )
        for col in extra.columns:
            meta[col] = extra[col]

        standard_base = base_cols if base_cols is not None else _STANDARD_OUTPUT_COLS
        extra_virus_cols = [c for c in extra.columns if c not in standard_base]
        out_cols = standard_base + extra_virus_cols
    else:
        meta.rename(columns=_CUSTOM_RENAME_MAP, inplace=True)
        if metadata_df is None:
            out_cols = ["Sequence_ID"]
        else:
            out_cols = _CUSTOM_OUTPUT_COLS

    for col in out_cols:
        if col not in meta.columns:
            meta[col] = ""

    sep = "\t" if not is_standard else ","
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if len(meta) <= NCBI_BATCH_SIZE:
        meta[out_cols].to_csv(out_path, sep=sep, index=False)
        logger.info("Metadata written: %s", out_path)
    else:
        for i in range(0, len(meta), NCBI_BATCH_SIZE):
            batch = meta.iloc[i : i + NCBI_BATCH_SIZE]
            batch_idx = i // NCBI_BATCH_SIZE + 1
            batch_path = _batch_path(out_path, batch_idx)
            batch[out_cols].to_csv(batch_path, sep=sep, index=False)
            logger.info(
                "Metadata batch %d written: %s (%d rows)",
                batch_idx,
                batch_path,
                len(batch),
            )
