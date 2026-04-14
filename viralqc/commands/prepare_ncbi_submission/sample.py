import logging
import pandas as pd
import re
import typer
from pathlib import Path
from typing import List, Optional

from viralqc.core.prepare_submission import PrepareSubmission

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


def _load_metadata_as_dicts(metadata_path: Path) -> list[dict]:
    """Read a metadata CSV and convert to the list-of-dicts format for PrepareSubmission."""
    df = pd.read_csv(metadata_path, dtype=str)
    col_remap = {
        "Sequence_ID": "sample_id",
        "geo_loc_name": "country",
        "host": "host",
        "isolate": "isolate",
        "collection-date": "collection-date",
        "isolation-source": "isolation-source",
    }
    df = df.rename(columns={k: v for k, v in col_remap.items() if k in df.columns})
    return df.to_dict(orient="records")


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
    metadata_dicts = _load_metadata_as_dicts(metadata)

    ps = PrepareSubmission(
        viralqc_results=results,
        viralqc_target_seq=sequences_vqc,
        viralqc_input_seq=sequences_input,
        samples_metadata=metadata_dicts,
        output_prefix=output_prefix,
        split_by_segments=split_by_segments,
        tbl_dir=tbl_dir,
    )

    entries = ps.run_sample(samples=flat_ids)

    total_sequences = sum(
        len(info["sequences"]) for e in entries for info in e.values()
    )

    summary_path = Path(f"{output_prefix}_summary.txt")
    skipped_path = Path(f"{output_prefix}_skipped.tsv")

    typer.echo(
        f"\nDone. {len(entries)} virus group(s) prepared, "
        f"{total_sequences} sequence file(s) total. "
        f"Details in {summary_path}."
        + (
            f" Skipped samples logged in {skipped_path}."
            if skipped_path.exists()
            else ""
        )
    )
