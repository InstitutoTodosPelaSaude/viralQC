import logging
import pandas as pd
import typer
from pathlib import Path
from typing import Optional
from viralqc.core.prepare_submission import PrepareSubmission

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


def _load_metadata_as_dicts(metadata_path: Path) -> list[dict]:
    """Read a metadata CSV and convert to list of dicts expected by PrepareSubmission."""
    df = pd.read_csv(metadata_path, dtype=str)
    # Columns in the CSV use the internal names (Sequence_ID, geo_loc_name, …).
    # Map them back to the user-facing key names so PrepareSubmission can re-map them.
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


def _make_submission(
    results: Path,
    sequences_vqc: Path,
    sequences_input: Path,
    output_prefix: str,
    metadata: Optional[Path],
    split_by_segments: bool = False,
    tbl_dir: Optional[Path] = None,
) -> PrepareSubmission:
    metadata_dicts = _load_metadata_as_dicts(metadata) if metadata else []
    return PrepareSubmission(
        viralqc_results=results,
        viralqc_target_seq=sequences_vqc,
        viralqc_input_seq=sequences_input,
        samples_metadata=metadata_dicts,
        output_prefix=output_prefix,
        split_by_segments=split_by_segments,
        tbl_dir=tbl_dir,
    )


@virus_app.command("sars-cov-2")
def sars_cov2(
    results: Path = _RESULTS_OPTION,
    sequences_vqc: Path = _SEQUENCES_VQC_OPTION,
    sequences_input: Path = _SEQUENCES_INPUT_OPTION,
    output_prefix: str = _PREFIX_OPTION,
    metadata: Path = _METADATA_OPTION,
):
    """Prepare NCBI submission package for SARS-CoV-2 sequences."""
    ps = _make_submission(
        results, sequences_vqc, sequences_input, output_prefix, metadata
    )
    entries = ps.run_virus(virus="sars-cov-2")
    for e in entries:
        for label, info in e.items():
            typer.echo(
                f"{label}: {len(info['sequences'])} file(s) — {info['sequences']}"
            )


@virus_app.command("dengue")
def dengue(
    results: Path = _RESULTS_OPTION,
    sequences_vqc: Path = _SEQUENCES_VQC_OPTION,
    sequences_input: Path = _SEQUENCES_INPUT_OPTION,
    output_prefix: str = _PREFIX_OPTION,
    metadata: Path = _METADATA_OPTION,
):
    """Prepare NCBI submission packages for Dengue virus (types 1–4)."""
    ps = _make_submission(
        results, sequences_vqc, sequences_input, output_prefix, metadata
    )
    entries = ps.run_virus(virus="dengue")
    for e in entries:
        for label, info in e.items():
            typer.echo(
                f"{label}: {len(info['sequences'])} file(s) — {info['sequences']}"
            )


@virus_app.command("influenza")
def influenza(
    results: Path = _RESULTS_OPTION,
    sequences_vqc: Path = _SEQUENCES_VQC_OPTION,
    sequences_input: Path = _SEQUENCES_INPUT_OPTION,
    output_prefix: str = _PREFIX_OPTION,
    metadata: Path = _METADATA_OPTION,
):
    """Prepare NCBI submission packages for Influenza A, B and C (per segment)."""
    ps = _make_submission(
        results, sequences_vqc, sequences_input, output_prefix, metadata
    )
    entries = ps.run_virus(virus="influenza")
    for e in entries:
        for label, info in e.items():
            typer.echo(
                f"{label}: {len(info['sequences'])} file(s) — {info['sequences']}"
            )


@virus_app.command("norovirus")
def norovirus(
    results: Path = _RESULTS_OPTION,
    sequences_vqc: Path = _SEQUENCES_VQC_OPTION,
    sequences_input: Path = _SEQUENCES_INPUT_OPTION,
    output_prefix: str = _PREFIX_OPTION,
    metadata: Path = _METADATA_OPTION,
):
    """Prepare NCBI submission packages for Norovirus (per genogroup GI–GVI)."""
    ps = _make_submission(
        results, sequences_vqc, sequences_input, output_prefix, metadata
    )
    entries = ps.run_virus(virus="norovirus")
    for e in entries:
        for label, info in e.items():
            typer.echo(
                f"{label}: {len(info['sequences'])} file(s) — {info['sequences']}"
            )


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
        "collection-date, isolation-source). Only Sequence_ID is required "
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
    split_by_segments: bool = typer.Option(
        False,
        "--split-by-segments",
        help="Organize output into per-segment subdirectories when the virus "
        "has annotated segments (e.g. S, M, L).",
    ),
):
    """Prepare an NCBI submission package for a custom (non-predefined) virus."""
    ps = _make_submission(
        results,
        sequences_vqc,
        sequences_input,
        output_prefix,
        metadata,
        split_by_segments=split_by_segments,
        tbl_dir=tbl_dir,
    )
    entries = ps.run_virus(virus="custom", virus_name=virus_name)
    for e in entries:
        for label, info in e.items():
            typer.echo(
                f"{label}: {len(info['sequences'])} file(s) — {info['sequences']}"
            )
