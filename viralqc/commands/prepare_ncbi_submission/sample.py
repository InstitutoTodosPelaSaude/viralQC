import logging
import pandas as pd
import re
import typer
from pathlib import Path
from typing import List, Optional
from viralqc.core.ncbi_submission import (
    filter_sequences,
    is_plain_header_virus,
    load_fasta,
    load_metadata,
    load_results,
    write_dropped_table,
    write_fasta,
    write_rename_log,
    write_submission_metadata,
    _NOROVIRUS_OUTPUT_COLS,
    _NOROVIRUS_RE,
)

logger = logging.getLogger(__name__)


def _collect_sample_ids(
    sample: Optional[List[str]],
    sample_ids: Optional[Path],
) -> List[str]:
    """
    Resolve the final list of sample IDs from CLI inputs.

    Args:
        sample: List of IDs passed via ``--sample`` (may contain 'all').
        sample_ids: Path to a text file with one ID per line (``--sample-ids``).

    Returns:
        Deduplicated list of sample ID strings, or the special value ``['all']``.
    """
    ids: List[str] = []

    if sample:
        for s in sample:
            ids.extend(re.split(r"[,\s]+", s.strip()))

    if sample_ids:
        text = sample_ids.read_text(encoding="utf-8")
        for line in text.splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                ids.append(line)

    ids = [s for s in ids if s]

    if "all" in ids:
        return ["all"]

    seen = set()
    return [x for x in ids if not (x in seen or seen.add(x))]


def _write_skipped_log(skipped: list, out_path: Path) -> None:
    """Write a TSV of samples that were skipped and why."""
    if not skipped:
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    import pandas as pd

    pd.DataFrame(skipped).to_csv(out_path, sep="\t", index=False)
    logger.warning("Some samples were skipped — see %s", out_path)


def _process_virus_group(
    virus: str,
    group_df,
    fasta: dict,
    tbl_dir: Optional[Path],
    output_prefix: str,
    skipped: list,
    summary_log: list,
    meta_df=None,
    split_by_segments: bool = False,
) -> None:
    """Build the submission package for a grouped virus slice."""
    from viralqc.commands.prepare_ncbi_submission.virus import _safe_dir_name

    sample_ids = group_df["seqName"].str.strip().tolist()

    if virus in ("", "nan", "Unclassified"):
        for sid in sample_ids:
            skipped.append(
                {
                    "sampleId": sid,
                    "reason": f"Unclassified virus ('{virus}')",
                }
            )
            logger.warning("Sample '%s' is unclassified — skipping.", sid)
        return 0

    use_plain_header = is_plain_header_virus(virus)

    if use_plain_header:
        if meta_df is None:
            typer.echo(
                f"Error: {virus} is a standard virus, but no --metadata CSV was provided. Metadata is required.",
                err=True,
            )
            raise typer.Exit(code=1)

        from viralqc.core.ncbi_submission import _STANDARD_OUTPUT_COLS

        req_cols = [
            c for c in _STANDARD_OUTPUT_COLS if c not in ("Sequence_ID", "serotype")
        ]
        missing = [c for c in req_cols if c not in meta_df.columns]
        if missing:
            typer.echo(
                f"Error: Provided metadata for {virus} is missing required columns: {missing}",
                err=True,
            )
            raise typer.Exit(code=1)

    group_seqs = {}
    for sid in sample_ids:
        found = False
        for k, v in fasta.items():
            if k == sid or k.startswith(sid):
                group_seqs[k] = v
                found = True
        if not found:
            skipped.append(
                {"sampleId": sid, "reason": "No sequence in target-regions FASTA"}
            )
            logger.warning(
                "Sample '%s' has no sequence in target-regions FASTA — skipping.", sid
            )

    if not group_seqs:
        return 0

    kept, dropped = filter_sequences(group_seqs, group_df)

    if not use_plain_header:
        valid_kept = []
        for record in kept:
            sid = record.id
            row = group_df[group_df["seqName"].str.strip() == sid].iloc[0]
            tbl_path_str = str(row.get("tbl_path", "") or "").strip()

            has_tbl = False
            if tbl_dir:
                candidates = list(tbl_dir.glob(f"*_{sid}.tbl")) or list(
                    tbl_dir.glob(f"*{re.escape(sid)}*.tbl")
                )
                if candidates and candidates[0].exists():
                    has_tbl = True
            else:
                if (
                    tbl_path_str
                    and tbl_path_str != "nan"
                    and Path(tbl_path_str).exists()
                ):
                    has_tbl = True

            if has_tbl:
                valid_kept.append(record)
            else:
                skipped.append(
                    {
                        "sampleId": sid,
                        "reason": f"TBL file missing or not found (virus: '{virus}')",
                    }
                )
        kept = valid_kept

    if not kept:
        return 0

    dir_virus_name = virus
    is_influenza = False
    is_norovirus = False
    if "SARS-CoV-2" in virus:
        dir_virus_name = "SARS-CoV-2"
    elif m := re.search(r"Dengue virus type\s+(\d)", virus, re.IGNORECASE):
        dir_virus_name = f"Dengue{m.group(1)}"
    elif m := re.search(r"^Influenza\s+([ABC])", virus, re.IGNORECASE):
        flu_type = m.group(1).upper()
        dir_virus_name = f"Influenza{flu_type}"
        is_influenza = True
    elif m := _NOROVIRUS_RE.search(virus):
        genogroup = m.group(1).upper()
        dir_virus_name = "Norovirus"
        is_norovirus = True
        noro_genogroup = genogroup

    base_out_dir = Path(f"{output_prefix}_{_safe_dir_name(dir_virus_name)}")
    base_out_dir.mkdir(parents=True, exist_ok=True)

    segment_groups = []
    if is_influenza:
        for seg_val, seg_df in group_df.groupby("segment", dropna=False):
            seg_str = str(seg_val).strip()
            seg_safe = (
                _safe_dir_name(seg_str) if seg_str and seg_str != "nan" else "unknown"
            )
            segment_groups.append((base_out_dir / seg_safe, seg_df))
    elif is_norovirus:
        segment_groups.append((base_out_dir / noro_genogroup, group_df))
    else:
        # Custom viruses: group by segment only when --split-by-segments is set
        has_segments = (
            split_by_segments
            and "segment" in group_df.columns
            and group_df["segment"].dropna().str.strip().replace("", pd.NA).dropna().any()
        )
        if has_segments:
            for seg_val, seg_df in group_df.groupby("segment", dropna=False):
                seg_str = str(seg_val).strip()
                seg_safe = (
                    _safe_dir_name(seg_str)
                    if seg_str and seg_str != "nan"
                    else "unknown"
                )
                segment_groups.append((base_out_dir / seg_safe, seg_df))
        else:
            segment_groups.append((base_out_dir, group_df))

    total_kept = 0
    for out_dir, seg_df in segment_groups:
        out_dir.mkdir(parents=True, exist_ok=True)
        seg_record_ids = set(seg_df["seqName"].str.strip())
        seg_kept = [r for r in kept if r.id in seg_record_ids]
        seg_dropped = [d for d in dropped if d["seqName"] in seg_record_ids]
        total_kept += len(seg_kept)

        rename_log: list = []
        virus_species = str(seg_df["virus_species"].iloc[0]).strip()

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
        )

        if meta_df is not None:
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

        if seg_dropped:
            write_dropped_table(seg_dropped, out_dir / "dropped_sequences.tsv")
            summary_log.append(
                f"{virus} [{out_dir.name}]: {len(seg_dropped)} sequence(s) dropped — see dropped_sequences.tsv"
            )

        if rename_log:
            write_rename_log(rename_log, out_dir / "renamed_headers.tsv")

        tbl_info = ""
        if not use_plain_header:
            tbl_parts: list[str] = []
            n_tbl = 0
            for record in seg_kept:
                sid = record.id
                row = seg_df[seg_df["seqName"].str.strip() == sid].iloc[0]
                tbl_content = None

                if tbl_dir:
                    candidates = list(tbl_dir.glob(f"*_{sid}.tbl")) or list(
                        tbl_dir.glob(f"*{re.escape(sid)}*.tbl")
                    )
                    if candidates and candidates[0].exists():
                        tbl_content = candidates[0].read_text(
                            encoding="utf-8", errors="replace"
                        )
                else:
                    tbl_path_str = str(row.get("tbl_path", "") or "").strip()
                    if (
                        tbl_path_str
                        and tbl_path_str != "nan"
                        and Path(tbl_path_str).exists()
                    ):
                        tbl_content = Path(tbl_path_str).read_text(
                            encoding="utf-8", errors="replace"
                        )

                if tbl_content:
                    tbl_parts.append(tbl_content.rstrip("\n"))
                    n_tbl += 1

            if tbl_parts:
                combined_tbl = out_dir / "annotation.tbl"
                combined_tbl.write_text(
                    "\n".join(tbl_parts) + "\n", encoding="utf-8"
                )
            if n_tbl > 0:
                tbl_info = f" + {n_tbl} TBL"

        summary_log.append(f"{virus}: {len(seg_kept)} seq(s){tbl_info} - {out_dir}")

    return total_kept


def sample_cmd(
    sample: Optional[List[str]] = typer.Option(
        None,
        "--sample",
        help=(
            "Sample ID to include. May be specified multiple times. "
            "Pass 'all' to process every sample in the results file."
        ),
    ),
    sample_ids: Optional[Path] = typer.Option(
        None,
        "--sample-ids",
        help="Path to a text file with one sample ID per line.",
    ),
    results: Path = typer.Option(
        ...,
        "--results",
        help="ViralQC results file (.tsv, .csv or .json).",
    ),
    sequences_vqc: Path = typer.Option(
        ...,
        "--sequences-vqc",
        help="ViralQC sequences_target_regions.fasta file.",
    ),
    sequences_input: Path = typer.Option(
        ...,
        "--sequences-input",
        help="Original input FASTA file provided to viralQC.",
    ),
    tbl_dir: Optional[Path] = typer.Option(
        None,
        "--tbl-dir",
        help="Directory containing per-sample TBL files. "
        "When omitted, uses the tbl_path column from the results file.",
    ),
    metadata: Path = typer.Option(
        ...,
        "--metadata",
        help="CSV with sample metadata (Sequence_ID, geo_loc_name, host, isolate, "
        "collection-date, isolation-source). When provided, a metadata file is "
        "written into each sample output directory.",
    ),
    output_prefix: str = typer.Option(
        "ncbi_submission",
        "--output-prefix",
        help="Prefix for the per-sample output directories.",
    ),
    split_by_segments: bool = typer.Option(
        False,
        "--split-by-segments",
        help="Organize custom virus output into per-segment subdirectories "
        "when the virus has annotated segments (e.g. S, M, L).",
    ),
):
    """Prepare NCBI submission packages for specific samples (or all samples)."""
    if not sample and not sample_ids:
        typer.echo("Provide at least --sample <ID> or --sample-ids <file>.", err=True)
        raise typer.Exit(code=1)

    flat_ids = _collect_sample_ids(sample, sample_ids)

    df = load_results(results)
    fasta = load_fasta(sequences_input)
    fasta.update(load_fasta(sequences_vqc))
    meta_df = load_metadata(metadata, is_standard=False)

    if flat_ids != ["all"]:
        df = df[df["seqName"].str.strip().isin(flat_ids)].copy()

    skipped: list = []
    summary_log: list = []
    total_kept = 0

    for virus_name, group_df in df.groupby("virus"):
        virus_str = str(virus_name).strip()
        n_kept = _process_virus_group(
            virus_str,
            group_df,
            fasta,
            tbl_dir,
            output_prefix,
            skipped,
            summary_log,
            meta_df,
            split_by_segments,
        )
        total_kept += n_kept or 0

    summary_path = Path(f"{output_prefix}_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("\n".join(summary_log) + "\n")

    skipped_unique = {(s["sampleId"], s["reason"]): s for s in skipped}
    _write_skipped_log(
        list(skipped_unique.values()), Path(f"{output_prefix}_skipped.tsv")
    )

    typer.echo(
        f"\nDone. {total_kept} sequences prepared for submission, "
        f"{len(skipped_unique)} sample(s) logged in {output_prefix}_skipped.tsv. "
        f"Detailed output summary written to {summary_path}"
    )
