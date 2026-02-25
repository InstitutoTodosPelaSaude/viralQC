import logging
import re
import typer
from pathlib import Path
from typing import Optional
from viralqc.core.ncbi_submission import (
    copy_tbl,
    filter_sequences,
    is_plain_header_virus,
    load_fasta,
    load_metadata,
    load_results,
    write_dropped_table,
    write_fasta,
    write_rename_log,
    write_submission_metadata,
    _SC2_OUTPUT_COLS,
)

logger = logging.getLogger(__name__)

virus_app = typer.Typer(
    name="virus",
    help="Prepare NCBI submission packages for a specific virus.",
    no_args_is_help=True,
)

_RESULTS_OPTION = typer.Option(
    ..., "--results", help="ViralQC results file (.tsv, .csv or .json)."
)
_SEQUENCES_VQC_OPTION = typer.Option(
    ..., "--sequences-vqc", help="ViralQC sequences_target_regions.fasta file."
)
_SEQUENCES_INPUT_OPTION = typer.Option(
    ..., "--sequences-input", help="Original input FASTA file provided to viralQC."
)
_PREFIX_OPTION = typer.Option(
    "ncbi_submission", "--output-prefix", help="Prefix for output directories."
)
_METADATA_OPTION = typer.Option(
    ...,
    "--metadata",
    help="CSV with sample metadata (Sequence_ID, geo_loc_name, host, isolate, "
    "collection-date, isolation-source).",
)


def common_virus_options(func):
    """Decorator that adds --results, --sequences, --output-prefix and --metadata to a command."""

    def wrapper(
        results: Path = _RESULTS_OPTION,
        sequences_vqc: Path = _SEQUENCES_VQC_OPTION,
        sequences_input: Path = _SEQUENCES_INPUT_OPTION,
        output_prefix: str = _PREFIX_OPTION,
        metadata: Path = _METADATA_OPTION,
    ):
        func(
            results=results,
            sequences_vqc=sequences_vqc,
            sequences_input=sequences_input,
            output_prefix=output_prefix,
            metadata=metadata,
        )

    return wrapper


def _safe_dir_name(s: str) -> str:
    """Replace whitespace and non-word chars with underscores for dir names."""
    return re.sub(r"[^\w.-]", "_", str(s)).strip("_")


def _build_package(
    out_dir: Path,
    records,
    dropped: list,
    rename_log: list,
    header_fn=None,
) -> None:
    """Write sequences.fasta, dropped_sequences.tsv and renamed_headers.tsv."""
    out_dir.mkdir(parents=True, exist_ok=True)

    write_fasta(
        records, out_dir / "sequences.fasta", header_fn=header_fn, rename_log=rename_log
    )

    if dropped:
        write_dropped_table(dropped, out_dir / "dropped_sequences.tsv")
        logger.warning(
            "%s: %d sequence(s) dropped — see dropped_sequences.tsv",
            out_dir.name,
            len(dropped),
        )

    if rename_log:
        write_rename_log(rename_log, out_dir / "renamed_headers.tsv")


@virus_app.command("sars-cov-2")
@common_virus_options
def sars_cov2(results, sequences_vqc, sequences_input, output_prefix, metadata):
    """Prepare NCBI submission package for SARS-CoV-2 sequences."""
    df = load_results(results)
    fasta = load_fasta(sequences_input)
    fasta.update(load_fasta(sequences_vqc))
    meta_df = load_metadata(metadata)

    sc2_mask = df["virus"].str.strip() == "SARS-CoV-2"
    sc2_seq_names = set(df.loc[sc2_mask, "seqName"].str.strip())

    subset = {k: v for k, v in fasta.items() if k in sc2_seq_names}
    kept, dropped = filter_sequences(subset, df)

    rename_log: list = []
    out_dir = Path(f"{output_prefix}_SARS-CoV-2")
    _build_package(out_dir, kept, dropped, rename_log, header_fn=lambda r: r.id)
    write_submission_metadata(
        kept,
        df,
        meta_df,
        out_dir / "metadata.csv",
        is_standard=True,
        base_cols=_SC2_OUTPUT_COLS,
    )

    typer.echo(f"SARS-CoV-2: {len(kept)} sequences - {out_dir}")


@virus_app.command("dengue")
@common_virus_options
def dengue(results, sequences_vqc, sequences_input, output_prefix, metadata):
    """Prepare NCBI submission packages for Dengue virus (types 1–4)."""
    df = load_results(results)
    fasta = load_fasta(sequences_input)
    fasta.update(load_fasta(sequences_vqc))
    meta_df = load_metadata(metadata)

    for dengue_type in (1, 2, 3, 4):
        virus_name = f"Dengue virus type {dengue_type}"
        mask = df["virus"].str.strip() == virus_name
        seq_names = set(df.loc[mask, "seqName"].str.strip())

        if not seq_names:
            logger.info("No sequences found for %s — skipping.", virus_name)
            continue

        subset = {k: v for k, v in fasta.items() if k in seq_names}
        kept, dropped = filter_sequences(subset, df)

        rename_log: list = []
        out_dir = Path(f"{output_prefix}_Dengue{dengue_type}")
        _build_package(out_dir, kept, dropped, rename_log, header_fn=lambda r: r.id)
        write_submission_metadata(
            kept, df, meta_df, out_dir / "metadata.csv", is_standard=True
        )

        typer.echo(f"Dengue type {dengue_type}: {len(kept)} sequences - {out_dir}")


def _influenza_group_key(row) -> tuple[str, str]:
    """
    Return a (type_dir, segment_dir) tuple for a given influenza result row.

    type_dir  — top-level directory name, e.g. ``'InfluenzaA_H1N1'``
    segment_dir — segment subdirectory name, e.g. ``'HA'``
    """
    virus = str(row.get("virus") or "")
    segment = str(row.get("segment") or "")
    seg_safe = _safe_dir_name(segment) if segment not in ("", "nan") else "unknown"

    m = re.match(r"Influenza\s+(A|B|C)", virus, re.IGNORECASE)
    if not m:
        return _safe_dir_name(virus), seg_safe

    flu_type = m.group(1).upper()
    return f"Influenza{flu_type}", seg_safe


@virus_app.command("influenza")
@common_virus_options
def influenza(results, sequences_vqc, sequences_input, output_prefix, metadata):
    """Prepare NCBI submission packages for Influenza A, B and C (per segment)."""
    df = load_results(results)
    fasta = load_fasta(sequences_input)
    fasta.update(load_fasta(sequences_vqc))
    meta_df = load_metadata(metadata)

    influenza_mask = df["virus"].str.contains(
        r"^Influenza\s+[ABC]", case=False, regex=True, na=False
    )
    df_flu = df[influenza_mask].copy()

    if df_flu.empty:
        typer.echo("No Influenza sequences found in results.")
        raise typer.Exit()

    df_flu["_type_dir"] = df_flu.apply(lambda r: _influenza_group_key(r)[0], axis=1)
    df_flu["_seg_dir"] = df_flu.apply(lambda r: _influenza_group_key(r)[1], axis=1)

    for (type_dir, seg_dir), group_df in df_flu.groupby(["_type_dir", "_seg_dir"]):
        seq_names = set(group_df["seqName"].str.strip())
        subset = {k: v for k, v in fasta.items() if k in seq_names}
        kept, dropped = filter_sequences(subset, df)

        rename_log: list = []
        out_dir = Path(output_prefix + "_" + type_dir) / seg_dir
        _build_package(out_dir, kept, dropped, rename_log, header_fn=lambda r: r.id)
        write_submission_metadata(
            kept, df, meta_df, out_dir / "metadata.csv", is_standard=True
        )

        typer.echo(f"{type_dir}/{seg_dir}: {len(kept)} sequences - {out_dir}")


@virus_app.command("custom")
def custom(
    results: Path = _RESULTS_OPTION,
    sequences_vqc: Path = _SEQUENCES_VQC_OPTION,
    sequences_input: Path = _SEQUENCES_INPUT_OPTION,
    output_prefix: str = _PREFIX_OPTION,
    metadata: Optional[Path] = typer.Option(
        None,
        "--metadata",
        help="Optional CSV with sample metadata (Sequence_ID, geo_loc_name, host, isolate, "
        "collection-date, isolation-source). Only Sequence_ID and Organism are required "
        "in the output; other columns are included when metadata is provided.",
    ),
    virus_name: str = typer.Option(
        ...,
        "--virus-name",
        help="Virus name used as the output directory suffix.",
    ),
    tbl_dir: Optional[Path] = typer.Option(
        None,
        "--tbl-dir",
        help="Directory containing per-sample TBL files. "
        "When omitted, uses the tbl_path column from the results file.",
    ),
):
    """Prepare an NCBI submission package for a custom (non-predefined) virus."""
    df = load_results(results)
    fasta = load_fasta(sequences_input)
    fasta.update(load_fasta(sequences_vqc))
    meta_df = load_metadata(metadata, is_standard=False) if metadata else None

    mask = ~df["virus"].apply(is_plain_header_virus)
    mask &= df["virus"].str.contains(re.escape(virus_name), case=False, na=False)
    df_custom = df[mask].copy()

    if df_custom.empty:
        typer.echo(f"No sequences found matching virus name '{virus_name}'.")
        raise typer.Exit()

    seq_names = set(df_custom["seqName"].str.strip())
    subset = {k: v for k, v in fasta.items() if k in seq_names}
    kept, dropped = filter_sequences(subset, df)

    row_lookup = {row["seqName"].strip(): row for _, row in df_custom.iterrows()}

    def header_fn(record):
        row = row_lookup.get(record.id, {})
        organism = row.get("virus_species", "")
        if not organism or str(organism) in ("nan", ""):
            organism = row.get("virus", virus_name)
        return f"{record.id} [Organism={organism}]"

    rename_log: list = []
    out_dir = Path(f"{output_prefix}_{_safe_dir_name(virus_name)}")
    _build_package(out_dir, kept, dropped, rename_log, header_fn=header_fn)

    first_row = row_lookup.get(kept[0].id, {}) if kept else {}
    organism_label = str(
        first_row.get("virus_species") or first_row.get("virus") or virus_name
    )
    write_submission_metadata(
        kept,
        df,
        meta_df,
        out_dir / "metadata.csv",
        is_standard=False,
        organism=organism_label,
    )

    n_tbl = 0
    tbl_out_dir = out_dir / "tbl_files"

    for record in kept:
        row = row_lookup.get(record.id, {})

        if tbl_dir:
            candidates = list(tbl_dir.glob(f"*_{record.id}.tbl"))
            if not candidates:
                candidates = list(tbl_dir.glob(f"*{re.escape(record.id)}*.tbl"))
            if candidates:
                if copy_tbl(
                    str(candidates[0]), tbl_out_dir, dest_name=f"{record.id}.tbl"
                ):
                    n_tbl += 1
        else:
            tbl_path_str = row.get("tbl_path", "")
            if copy_tbl(str(tbl_path_str), tbl_out_dir, dest_name=f"{record.id}.tbl"):
                n_tbl += 1

    typer.echo(f"{virus_name}: {len(kept)} sequences, {n_tbl} TBL files - {out_dir}")
