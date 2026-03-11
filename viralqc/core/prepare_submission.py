from __future__ import annotations
import logging
import re
import pandas as pd
from pathlib import Path
from typing import List, Optional
from viralqc.core.ncbi_submission import (
    _NOROVIRUS_OUTPUT_COLS,
    _NOROVIRUS_RE,
    _SC2_OUTPUT_COLS,
    _STANDARD_OUTPUT_COLS,
    filter_sequences,
    is_plain_header_virus,
    load_fasta,
    load_results,
    write_fasta,
    write_rename_log,
    write_submission_metadata,
)

logger = logging.getLogger(__name__)

_METADATA_KEY_MAP = {
    "sample_id": "Sequence_ID",
    "country": "geo_loc_name",
    "host": "host",
    "isolate": "isolate",
    "collection-date": "collection-date",
    "isolation-source": "isolation-source",
}


def _safe_dir_name(s: str) -> str:
    """Replace whitespace and non-word chars with underscores for dir names."""
    return re.sub(r"[^\w.-]", "_", str(s)).strip("_")


def _norovirus_genogroup(virus: str) -> str:
    m = _NOROVIRUS_RE.search(virus)
    return m.group(1).upper() if m else "unknown"


def _influenza_group_key(row: dict) -> tuple[str, str]:
    virus = str(row.get("virus") or "")
    segment = str(row.get("segment") or "")
    seg_safe = _safe_dir_name(segment) if segment not in ("", "nan") else "unknown"
    m = re.match(r"Influenza\s+(A|B|C)", virus, re.IGNORECASE)
    if not m:
        return _safe_dir_name(virus), seg_safe
    flu_type = m.group(1).upper()
    return f"Influenza{flu_type}", seg_safe


def _collect_fasta_paths_from_dir(out_dir: Path) -> list[Path]:
    return sorted(out_dir.glob("sequences*.fasta"))


def _collect_metadata_paths_from_dir(out_dir: Path) -> list[Path]:
    return sorted(out_dir.glob("metadata*.*"))


def _collect_annotation_paths_from_dir(out_dir: Path) -> list[Path]:
    tbls = sorted(out_dir.glob("annotation*.tbl"))
    return tbls if tbls else []


def _build_package(
    out_dir: Path,
    records: list,
    dropped: list,
    rename_log: list,
    header_fn=None,
) -> None:
    """Write sequences.fasta and renamed_headers.tsv."""
    out_dir.mkdir(parents=True, exist_ok=True)
    write_fasta(
        records,
        out_dir / "sequences.fasta",
        header_fn=header_fn,
        rename_log=rename_log,
        dropped_log=dropped,
    )
    if rename_log:
        write_rename_log(rename_log, out_dir / "renamed_headers.tsv")


def _entry_for_dir(virus_label: str, out_dir: Path) -> dict:
    """Build one result entry dict for a given output directory."""
    entry: dict = {
        "sequences": _collect_fasta_paths_from_dir(out_dir),
        "metadata": _collect_metadata_paths_from_dir(out_dir),
        "log": out_dir / "summary.txt",
    }
    annotations = _collect_annotation_paths_from_dir(out_dir)
    if annotations:
        entry["annotation"] = annotations
    return {virus_label: entry}


class PrepareSubmission:
    """Prepare NCBI submission packages programmatically.

    This class mirrors the ``vqc prepare-ncbi-submission`` CLI behaviour but
    accepts metadata as a list of dicts and returns file paths instead of
    printing to stdout.

    Args:
        viralqc_results: Path to the ViralQC results file (.tsv, .csv or .json).
        viralqc_target_seq: Path to the target-regions FASTA produced by viralQC run.
        viralqc_input_seq: Path to the original input FASTA passed to viralQC run.
        samples_metadata: List of sample metadata dicts. Each dict must contain
            ``sample_id`` (≤ 24 chars). Standard viruses additionally require
            ``country``, ``host``, ``isolate``, ``collection-date`` and
            ``isolation-source``.
        output_prefix: Prefix for all output directories (default: ``ncbi_submission``).
        split_by_segments: When ``True``, custom viruses with annotated segments are
            split into per-segment subdirectories.
        tbl_dir: Optional directory containing per-sample ``.tbl`` annotation files.
            When *None*, the ``tbl_path`` column from the results file is used.
    """

    def __init__(
        self,
        viralqc_results: Path,
        viralqc_target_seq: Path,
        viralqc_input_seq: Path,
        samples_metadata: List[dict],
        output_prefix: str = "ncbi_submission",
        split_by_segments: bool = False,
        tbl_dir: Optional[Path] = None,
    ) -> None:
        self.viralqc_results = Path(viralqc_results)
        self.viralqc_target_seq = Path(viralqc_target_seq)
        self.viralqc_input_seq = Path(viralqc_input_seq)
        self.samples_metadata = samples_metadata
        self.output_prefix = output_prefix
        self.split_by_segments = split_by_segments
        self.tbl_dir = Path(tbl_dir) if tbl_dir else None

    def run_virus(
        self,
        virus: str = "all",
        virus_name: Optional[str] = None,
    ) -> List[dict]:
        """Prepare submission packages grouped by virus type.

        Mirrors ``vqc prepare-ncbi-submission virus <subcommand>``.

        When *virus* is ``"all"`` (default), every virus found in the results
        file is processed, applying the correct logic for each (standard or
        custom). You can also pass one of the predefined values: ``"sars-cov-2"``,
        ``"dengue"``, ``"influenza"``, ``"norovirus"``, or ``"custom"`` (requires
        *virus_name*).

        Args:
            virus: Which virus group to process. ``"all"`` processes everything.
            virus_name: Required when *virus* is ``"custom"``; the name to match
                against the ``virus`` column of the results file.

        Returns:
            A list of dicts, one per generated directory::

                [
                    {"SARS-CoV-2": {"sequences": [...], "metadata": [...], "log": Path, "annotation": [...]}},
                    {"Dengue1":    {"sequences": [...], ...}},
                ]
        """
        df = load_results(self.viralqc_results)
        fasta = load_fasta(self.viralqc_input_seq)
        fasta.update(load_fasta(self.viralqc_target_seq))
        meta_df = self._metadata_df()

        results: list[dict] = []

        if virus == "custom":
            if not virus_name:
                raise ValueError("'virus_name' is required when virus='custom'.")
            results.extend(self._process_custom_virus(df, fasta, meta_df, virus_name))
            return results

        if virus in ("sars-cov-2", "all"):
            results.extend(self._process_sars_cov2(df, fasta, meta_df))

        if virus in ("dengue", "all"):
            results.extend(self._process_dengue(df, fasta, meta_df))

        if virus in ("influenza", "all"):
            results.extend(self._process_influenza(df, fasta, meta_df))

        if virus in ("norovirus", "all"):
            results.extend(self._process_norovirus(df, fasta, meta_df))

        if virus == "all":
            # Custom viruses: everything that isn't a plain-header (standard) virus
            results.extend(self._process_all_custom(df, fasta, meta_df))

        return results

    def run_sample(
        self,
        samples: List[str] = ["all"],
    ) -> List[dict]:
        """Prepare submission packages for specific samples (or all).

        Mirrors ``vqc prepare-ncbi-submission sample``.

        Args:
            samples: List of sample IDs to process. Pass ``["all"]`` to process
                every sample present in the results file.

        Returns:
            Same list-of-dicts structure as :meth:`run_virus`.
        """
        df = load_results(self.viralqc_results)
        fasta = load_fasta(self.viralqc_input_seq)
        fasta.update(load_fasta(self.viralqc_target_seq))
        meta_df = self._metadata_df()

        if samples != ["all"]:
            flat_ids = set(samples)
            df = df[df["seqName"].str.strip().isin(flat_ids)].copy()

        results: list[dict] = []
        skipped: list[dict] = []
        summary_lines: list[str] = []

        for virus_name, group_df in df.groupby("virus"):
            virus_str = str(virus_name).strip()
            entries = self._process_sample_virus_group(
                virus_str, group_df, fasta, meta_df, skipped, summary_lines
            )
            results.extend(entries)

        # Write summary log
        summary_path = Path(f"{self.output_prefix}_summary.txt")
        with open(summary_path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(summary_lines) + "\n")

        # Write skipped log
        if skipped:
            skipped_unique = {(s["sampleId"], s["reason"]): s for s in skipped}
            skipped_path = Path(f"{self.output_prefix}_skipped.tsv")
            pd.DataFrame(list(skipped_unique.values())).to_csv(
                skipped_path, sep="\t", index=False
            )
            logger.warning("Some samples skipped — see %s", skipped_path)

        return results

    def _metadata_df(self) -> pd.DataFrame:
        """Convert samples_metadata list of dicts to a validated DataFrame."""
        if not self.samples_metadata:
            return pd.DataFrame(columns=["Sequence_ID"])

        rows = []
        for item in self.samples_metadata:
            row = {_METADATA_KEY_MAP.get(k, k): v for k, v in item.items()}
            rows.append(row)

        df = pd.DataFrame(rows, dtype=str)

        if "Sequence_ID" not in df.columns:
            raise ValueError(
                "Each metadata dict must contain a 'sample_id' key (maps to Sequence_ID)."
            )

        long_ids = df[df["Sequence_ID"].str.len() >= 25]["Sequence_ID"].tolist()
        if long_ids:
            raise ValueError(
                f"Sequence_ID values must be less than 25 characters. "
                f"Found {len(long_ids)} invalid IDs: {long_ids[:5]}"
            )

        return df

    def _process_sars_cov2(self, df, fasta, meta_df) -> list[dict]:
        mask = df["virus"].str.strip() == "SARS-CoV-2"
        seq_names = set(df.loc[mask, "seqName"].str.strip())
        if not seq_names:
            return []

        subset = {k: v for k, v in fasta.items() if k in seq_names}
        kept, dropped = filter_sequences(subset, df)
        if not kept:
            return []

        out_dir = Path(f"{self.output_prefix}_SARS-CoV-2")
        rename_log: list = []
        _build_package(out_dir, kept, dropped, rename_log, header_fn=lambda r: r.id)
        write_submission_metadata(
            kept,
            df,
            meta_df,
            out_dir / "metadata.csv",
            is_standard=True,
            base_cols=_SC2_OUTPUT_COLS,
        )
        self._write_summary(out_dir, f"SARS-CoV-2: {len(kept)} sequences", dropped)
        return [_entry_for_dir("SARS-CoV-2", out_dir)]

    def _process_dengue(self, df, fasta, meta_df) -> list[dict]:
        entries: list[dict] = []
        for dengue_type in (1, 2, 3, 4):
            virus_name = f"Dengue virus type {dengue_type}"
            mask = df["virus"].str.strip() == virus_name
            seq_names = set(df.loc[mask, "seqName"].str.strip())
            if not seq_names:
                continue

            subset = {k: v for k, v in fasta.items() if k in seq_names}
            kept, dropped = filter_sequences(subset, df)
            if not kept:
                continue

            out_dir = Path(f"{self.output_prefix}_Dengue{dengue_type}")
            rename_log: list = []
            _build_package(out_dir, kept, dropped, rename_log, header_fn=lambda r: r.id)
            write_submission_metadata(
                kept, df, meta_df, out_dir / "metadata.csv", is_standard=True
            )
            self._write_summary(
                out_dir, f"Dengue{dengue_type}: {len(kept)} sequences", dropped
            )
            entries.append(_entry_for_dir(f"Dengue{dengue_type}", out_dir))
        return entries

    def _process_influenza(self, df, fasta, meta_df) -> list[dict]:
        influenza_mask = df["virus"].str.contains(
            r"^Influenza\s+[ABC]", case=False, regex=True, na=False
        )
        df_flu = df[influenza_mask].copy()
        if df_flu.empty:
            return []

        df_flu["_type_dir"] = df_flu.apply(lambda r: _influenza_group_key(r)[0], axis=1)
        df_flu["_seg_dir"] = df_flu.apply(lambda r: _influenza_group_key(r)[1], axis=1)

        entries: list[dict] = []
        for (type_dir, seg_dir), group_df in df_flu.groupby(["_type_dir", "_seg_dir"]):
            seq_names = set(group_df["seqName"].str.strip())
            subset = {k: v for k, v in fasta.items() if k in seq_names}
            kept, dropped = filter_sequences(subset, df)
            if not kept:
                continue

            out_dir = Path(self.output_prefix + "_" + type_dir) / seg_dir
            rename_log: list = []
            _build_package(out_dir, kept, dropped, rename_log, header_fn=lambda r: r.id)
            write_submission_metadata(
                kept, df, meta_df, out_dir / "metadata.csv", is_standard=True
            )
            label = f"{type_dir}/{seg_dir}"
            self._write_summary(out_dir, f"{label}: {len(kept)} sequences", dropped)
            entries.append(_entry_for_dir(label, out_dir))
        return entries

    def _process_norovirus(self, df, fasta, meta_df) -> list[dict]:
        noro_mask = df["virus"].str.contains(
            r"^Norovirus", case=False, regex=True, na=False
        )
        df_noro = df[noro_mask].copy()
        if df_noro.empty:
            return []

        df_noro["_genogroup"] = df_noro["virus"].apply(_norovirus_genogroup)

        entries: list[dict] = []
        for genogroup, group_df in df_noro.groupby("_genogroup"):
            seq_names = set(group_df["seqName"].str.strip())
            subset = {k: v for k, v in fasta.items() if k in seq_names}
            kept, dropped = filter_sequences(subset, df)
            if not kept:
                continue

            out_dir = Path(f"{self.output_prefix}_Norovirus") / genogroup
            rename_log: list = []
            _build_package(out_dir, kept, dropped, rename_log, header_fn=lambda r: r.id)
            write_submission_metadata(
                kept,
                df,
                meta_df,
                out_dir / "metadata.csv",
                is_standard=True,
                base_cols=_NOROVIRUS_OUTPUT_COLS,
            )
            label = f"Norovirus/{genogroup}"
            self._write_summary(out_dir, f"{label}: {len(kept)} sequences", dropped)
            entries.append(_entry_for_dir(label, out_dir))
        return entries

    def _process_custom_virus(
        self, df, fasta, meta_df: Optional[pd.DataFrame], virus_name: str
    ) -> list[dict]:
        mask = ~df["virus"].apply(is_plain_header_virus)
        mask &= df["virus"].str.contains(re.escape(virus_name), case=False, na=False)
        df_custom = df[mask].copy()

        if df_custom.empty:
            logger.info("No sequences found matching virus name '%s'.", virus_name)
            return []

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

        base_out_dir = Path(f"{self.output_prefix}_{_safe_dir_name(virus_name)}")

        has_segments = (
            self.split_by_segments
            and "segment" in df_custom.columns
            and df_custom["segment"]
            .dropna()
            .str.strip()
            .replace("", pd.NA)
            .dropna()
            .any()
        )

        if has_segments:
            segment_groups = []
            for seg_val, seg_df in df_custom.groupby("segment", dropna=False):
                seg_str = str(seg_val).strip()
                seg_safe = (
                    _safe_dir_name(seg_str)
                    if seg_str and seg_str != "nan"
                    else "unknown"
                )
                seg_ids = set(seg_df["seqName"].str.strip())
                seg_kept = [r for r in kept if r.id in seg_ids]
                seg_dropped = [d for d in dropped if d["seqName"] in seg_ids]
                segment_groups.append(
                    (base_out_dir / seg_safe, seg_kept, seg_dropped, seg_df)
                )
        else:
            segment_groups = [(base_out_dir, kept, dropped, df_custom)]

        entries: list[dict] = []
        for out_dir, seg_kept, seg_dropped, seg_df in segment_groups:
            if not seg_kept:
                continue

            rename_log: list = []
            _build_package(
                out_dir, seg_kept, seg_dropped, rename_log, header_fn=header_fn
            )

            first_row = row_lookup.get(seg_kept[0].id, {}) if seg_kept else {}
            organism_label = str(
                first_row.get("virus_species") or first_row.get("virus") or virus_name
            )

            if meta_df is not None and not meta_df.empty:
                write_submission_metadata(
                    seg_kept,
                    df,
                    meta_df,
                    out_dir / "metadata.tsv",
                    is_standard=False,
                    organism=organism_label,
                )

            tbl_parts = self._collect_tbl(seg_kept, seg_df, row_lookup)
            if tbl_parts:
                (out_dir / "annotation.tbl").write_text(
                    "\n".join(tbl_parts) + "\n", encoding="utf-8"
                )

            label = f"{virus_name}/{out_dir.name}" if has_segments else virus_name
            self._write_summary(
                out_dir, f"{label}: {len(seg_kept)} sequences", seg_dropped
            )
            entries.append(_entry_for_dir(label, out_dir))

        return entries

    def _process_all_custom(self, df, fasta, meta_df) -> list[dict]:
        """Process all non-standard viruses found in results."""
        custom_mask = ~df["virus"].apply(is_plain_header_virus)
        df_custom = df[custom_mask].copy()
        if df_custom.empty:
            return []

        entries: list[dict] = []
        for virus_name, group_df in df_custom.groupby("virus"):
            vname = str(virus_name).strip()
            if vname in ("", "nan", "Unclassified"):
                continue
            entries.extend(self._process_custom_virus(df, fasta, meta_df, vname))
        return entries

    def _process_sample_virus_group(
        self,
        virus: str,
        group_df,
        fasta: dict,
        meta_df: pd.DataFrame,
        skipped: list,
        summary_lines: list,
    ) -> list[dict]:
        """Process one virus slice during run_sample(); mirrors sample.py logic."""
        if virus in ("", "nan", "Unclassified"):
            for sid in group_df["seqName"].str.strip().tolist():
                skipped.append(
                    {"sampleId": sid, "reason": f"Unclassified virus ('{virus}')"}
                )
            return []

        use_plain_header = is_plain_header_virus(virus)

        # Validate that standard-virus metadata is present
        if use_plain_header and meta_df is not None:
            req_cols = [
                c for c in _STANDARD_OUTPUT_COLS if c not in ("Sequence_ID", "serotype")
            ]
            missing = [c for c in req_cols if c not in meta_df.columns]
            if missing:
                logger.error(
                    "Metadata for %s is missing required columns: %s — skipping.",
                    virus,
                    missing,
                )
                return []

        sample_ids = group_df["seqName"].str.strip().tolist()
        group_seqs: dict = {}
        for sid in sample_ids:
            for k, v in fasta.items():
                if k == sid or k.startswith(sid):
                    group_seqs[k] = v
                    break
            else:
                skipped.append(
                    {"sampleId": sid, "reason": "No sequence in target-regions FASTA"}
                )

        if not group_seqs:
            return []

        kept, dropped = filter_sequences(group_seqs, group_df)

        # For custom viruses, skip sequences without TBL
        if not use_plain_header:
            valid_kept = []
            for record in kept:
                if self._has_tbl(record.id, group_df):
                    valid_kept.append(record)
                else:
                    skipped.append(
                        {
                            "sampleId": record.id,
                            "reason": f"TBL file missing or not found (virus: '{virus}')",
                        }
                    )
            kept = valid_kept

        if not kept:
            return []

        # Determine output directory name and grouping
        dir_virus_name, is_influenza, is_norovirus, noro_genogroup = (
            self._classify_virus(virus)
        )
        base_out_dir = Path(f"{self.output_prefix}_{_safe_dir_name(dir_virus_name)}")

        segment_groups = self._build_segment_groups(
            base_out_dir,
            group_df,
            kept,
            dropped,
            virus,
            is_influenza,
            is_norovirus,
            noro_genogroup if is_norovirus else None,
        )

        row_lookup = {row["seqName"].strip(): row for _, row in group_df.iterrows()}
        entries: list[dict] = []

        for out_dir, seg_df, seg_kept, seg_dropped in segment_groups:
            out_dir.mkdir(parents=True, exist_ok=True)
            rename_log: list = []
            virus_species = (
                str(seg_df["virus_species"].iloc[0]).strip()
                if "virus_species" in seg_df.columns
                else ""
            )

            if use_plain_header:

                def header_fn(r):
                    return r.id

            else:
                organism = virus_species if virus_species not in ("nan", "") else virus

                def header_fn(r, _org=organism):
                    return f"{r.id} [Organism={_org}]"

            write_fasta(
                seg_kept,
                out_dir / "sequences.fasta",
                header_fn=header_fn,
                rename_log=rename_log,
                dropped_log=seg_dropped,
            )

            if meta_df is not None and not meta_df.empty:
                organism_label = (
                    virus_species if virus_species not in ("nan", "") else virus
                )
                meta_ext = "csv" if use_plain_header else "tsv"
                noro_base = _NOROVIRUS_OUTPUT_COLS if is_norovirus else None
                write_submission_metadata(
                    seg_kept,
                    seg_df,
                    meta_df,
                    out_dir / f"metadata.{meta_ext}",
                    is_standard=use_plain_header,
                    organism=organism_label,
                    base_cols=noro_base,
                )

            if rename_log:
                write_rename_log(rename_log, out_dir / "renamed_headers.tsv")

            if not use_plain_header:
                tbl_parts = self._collect_tbl(seg_kept, seg_df, row_lookup)
                if tbl_parts:
                    (out_dir / "annotation.tbl").write_text(
                        "\n".join(tbl_parts) + "\n", encoding="utf-8"
                    )

            label = f"{virus}/{out_dir.name}" if len(segment_groups) > 1 else virus
            line = f"{label}: {len(seg_kept)} seq(s)"
            summary_lines.append(line)
            self._write_summary(out_dir, line, seg_dropped)
            entries.append(_entry_for_dir(label, out_dir))

        return entries

    def _classify_virus(self, virus: str) -> tuple[str, bool, bool, str]:
        """Return (dir_name, is_influenza, is_norovirus, noro_genogroup)."""
        is_influenza = False
        is_norovirus = False
        noro_genogroup = ""
        dir_name = virus

        if "SARS-CoV-2" in virus:
            dir_name = "SARS-CoV-2"
        elif m := re.search(r"Dengue virus type\s+(\d)", virus, re.IGNORECASE):
            dir_name = f"Dengue{m.group(1)}"
        elif re.search(r"^Influenza\s+[ABC]", virus, re.IGNORECASE):
            m2 = re.search(r"^Influenza\s+([ABC])", virus, re.IGNORECASE)
            flu_type = m2.group(1).upper() if m2 else ""
            dir_name = f"Influenza{flu_type}"
            is_influenza = True
        elif m3 := _NOROVIRUS_RE.search(virus):
            dir_name = "Norovirus"
            is_norovirus = True
            noro_genogroup = m3.group(1).upper()

        return dir_name, is_influenza, is_norovirus, noro_genogroup

    def _build_segment_groups(
        self,
        base_out_dir: Path,
        group_df,
        kept: list,
        dropped: list,
        virus: str,
        is_influenza: bool,
        is_norovirus: bool,
        noro_genogroup: Optional[str],
    ) -> list[tuple]:
        """Return list of (out_dir, seg_df, seg_kept, seg_dropped)."""
        groups = []
        if is_influenza:
            for seg_val, seg_df in group_df.groupby("segment", dropna=False):
                seg_str = str(seg_val).strip()
                seg_safe = (
                    _safe_dir_name(seg_str)
                    if seg_str and seg_str != "nan"
                    else "unknown"
                )
                out_dir = base_out_dir / seg_safe
                seg_ids = set(seg_df["seqName"].str.strip())
                seg_kept = [r for r in kept if r.id in seg_ids]
                seg_dropped = [d for d in dropped if d["seqName"] in seg_ids]
                groups.append((out_dir, seg_df, seg_kept, seg_dropped))
        elif is_norovirus:
            out_dir = base_out_dir / (noro_genogroup or "unknown")
            groups.append((out_dir, group_df, kept, dropped))
        else:
            has_segs = (
                self.split_by_segments
                and "segment" in group_df.columns
                and group_df["segment"]
                .dropna()
                .str.strip()
                .replace("", pd.NA)
                .dropna()
                .any()
            )
            if has_segs:
                for seg_val, seg_df in group_df.groupby("segment", dropna=False):
                    seg_str = str(seg_val).strip()
                    seg_safe = (
                        _safe_dir_name(seg_str)
                        if seg_str and seg_str != "nan"
                        else "unknown"
                    )
                    out_dir = base_out_dir / seg_safe
                    seg_ids = set(seg_df["seqName"].str.strip())
                    seg_kept = [r for r in kept if r.id in seg_ids]
                    seg_dropped = [d for d in dropped if d["seqName"] in seg_ids]
                    groups.append((out_dir, seg_df, seg_kept, seg_dropped))
            else:
                groups.append((base_out_dir, group_df, kept, dropped))
        return groups

    def _has_tbl(self, sid: str, group_df) -> bool:
        """Return True if a TBL file exists for the given sample ID."""
        if self.tbl_dir:
            candidates = list(self.tbl_dir.glob(f"*_{sid}.tbl")) or list(
                self.tbl_dir.glob(f"*{re.escape(sid)}*.tbl")
            )
            return bool(candidates and candidates[0].exists())
        # Read tbl_path from the results row
        rows = group_df[group_df["seqName"].str.strip() == sid]
        if rows.empty:
            return False
        tbl_path_str = str(rows.iloc[0].get("tbl_path", "") or "").strip()
        return bool(
            tbl_path_str and tbl_path_str != "nan" and Path(tbl_path_str).exists()
        )

    def _collect_tbl(self, seg_kept: list, seg_df, row_lookup: dict) -> list[str]:
        """Collect TBL file contents for all records in seg_kept."""
        parts: list[str] = []
        for record in seg_kept:
            sid = record.id
            content = None
            if self.tbl_dir:
                candidates = list(self.tbl_dir.glob(f"*_{sid}.tbl")) or list(
                    self.tbl_dir.glob(f"*{re.escape(sid)}*.tbl")
                )
                if candidates and candidates[0].exists():
                    content = candidates[0].read_text(
                        encoding="utf-8", errors="replace"
                    )
            else:
                row = row_lookup.get(sid, {})
                tbl_path_str = str(row.get("tbl_path", "") or "").strip()
                if (
                    tbl_path_str
                    and tbl_path_str != "nan"
                    and Path(tbl_path_str).exists()
                ):
                    content = Path(tbl_path_str).read_text(
                        encoding="utf-8", errors="replace"
                    )
            if content:
                parts.append(content.rstrip("\n"))
        return parts

    @staticmethod
    def _write_summary(
        out_dir: Path, message: str, dropped: list | None = None
    ) -> None:
        """Write summary.txt with the main message and, if any, dropped-sequence details."""
        out_dir.mkdir(parents=True, exist_ok=True)
        lines = [message]
        if dropped:
            lines.append(f"  dropped: {len(dropped)} sequence(s)")
            for item in dropped:
                seq = item.get("seqName", "?")
                reason = item.get("reason", "")
                n_pct = item.get("n_pct", "")
                length = item.get("length", "")
                detail = f"    - {seq}: {reason}"
                if length:
                    detail += f" (length={length}"
                    if n_pct:
                        detail += f", N%={n_pct}"
                    detail += ")"
                lines.append(detail)
        (out_dir / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
