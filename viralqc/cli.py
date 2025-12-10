import typer, logging, colorlog
from typing import Optional
from enum import Enum
from viralqc.core.datasets import GetNextcladeDatasets, GetBlastDatabase
from viralqc.core.run_analysis import RunAnalysis
from viralqc.core.defaults import (
    DEFAULT_SNK_FILE,
    DEFAULT_CONFIG_FILE,
    DEFAULT_CORES,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_OUTPUT_FILE,
    DEFAULT_DATASETS_LOCAL_PATH,
    DEFAULT_NEXTCLADE_SORT_MIN_SCORE,
    DEFAULT_NEXTCLADE_SORT_MIN_HITS,
    DEFAULT_BLAST_DATABASE,
    DEFAULT_BLAST_DATABASE_METADATA,
    DEFAULT_BLAST_IDENTITY_THRESHOLD,
    DEFAULT_BLAST_EVALUE,
    DEFAULT_BLAST_QCOV,
    DEFAULT_BLAST_TASK,
    DEFAULT_VERBOSE,
)
from viralqc import (
    GET_NC_PUBLIC_DATASETS_SNK_PATH,
    GET_BLAST_DB_SNK_PATH,
)

# core config
get_nc_datasets = GetNextcladeDatasets()
get_blast_db = GetBlastDatabase()
run_analysis = RunAnalysis()

# log config
handler = colorlog.StreamHandler()
handler.setLevel(logging.DEBUG)
handler.setFormatter(
    colorlog.ColoredFormatter(
        "%(log_color)s%(levelname)s:%(name)s:%(message)s",
        log_colors={
            "INFO": "green",
            "WARNING": "yellow",
            "ERROR": "red",
        },
    )
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(handler)


def log_multiline(text: str):
    for line in text.splitlines():
        if line.startswith("WARNING:"):
            logger.warning(line[len("WARNING:") :].strip())
        elif line.startswith("ERROR:"):
            logger.error(line[len("ERROR:") :].strip())
        elif line.startswith("INFO:"):
            logger.info(line[len("INFO:") :].strip())
        else:
            logger.debug(line.strip())


# cli config
app = typer.Typer()
if __name__ == "__main__":
    app()


@app.command()
def get_nextclade_datasets(
    datasets_dir: str = typer.Option(
        "datasets",
        "--datasets-dir",
        help="Directory to store local nextclade datasets.",
    ),
    snk_file_path: Optional[str] = GET_NC_PUBLIC_DATASETS_SNK_PATH,
    config_file_path: Optional[str] = DEFAULT_CONFIG_FILE,
    cores: int = 1,
    verbose: bool = typer.Option(
        DEFAULT_VERBOSE, "--verbose", "-v", help="Show snakemake logs."
    ),
):
    """Get Nextclade virus datasets"""
    snakemake_response = get_nc_datasets.get_public_dataset(
        datasets_dir=datasets_dir,
        snk_file=snk_file_path,
        config_file=config_file_path,
        cores=cores,
        verbose=verbose,
    )
    if snakemake_response.status == 200:
        if verbose:
            log_multiline(snakemake_response.format_log())
        logger.info("Nextclade public datasets successfully retrieved.")
    else:
        log_multiline(snakemake_response.format_log())
        logger.error("Failed to retrieve Nextclade public datasets.")


def validate_date_format(value: Optional[str]) -> Optional[str]:
    """Validate that the date is in YYYY-MM-DD format."""
    if value is None:
        return None
    import re
    from datetime import datetime

    if not re.match(r"^\d{4}-\d{2}-\d{2}$", value):
        raise typer.BadParameter("Date must be in YYYY-MM-DD format.")
    try:
        datetime.strptime(value, "%Y-%m-%d")
    except ValueError:
        raise typer.BadParameter(
            "Invalid date. Please use a valid date in YYYY-MM-DD format."
        )
    return value


@app.command()
def get_blast_database(
    output_dir: str = typer.Option(
        "datasets",
        "--output-dir",
        help="Path to store BLAST database.",
    ),
    release_date: Optional[str] = typer.Option(
        None,
        "--release-date",
        help="Filter sequences by release date (YYYY-MM-DD). Only sequences released on or before this date will be included.",
        callback=validate_date_format,
    ),
    snk_file_path: Optional[str] = GET_BLAST_DB_SNK_PATH,
    cores: int = 1,
    verbose: bool = typer.Option(
        DEFAULT_VERBOSE, "--verbose", "-v", help="Show snakemake logs."
    ),
):
    """Create BLAST database based on ncbi viruses refseq genomes"""
    snakemake_response = get_blast_db.get_database(
        output_dir=output_dir,
        release_date=release_date,
        snk_file=snk_file_path,
        cores=cores,
        verbose=verbose,
    )
    if snakemake_response.status == 200:
        if verbose:
            log_multiline(snakemake_response.format_log())
        logger.info("BLAST database created.")
    else:
        log_multiline(snakemake_response.format_log())
        logger.error("Failed to create BLAST database.")


class SortChoices(str, Enum):
    nextclade = ("nextclade",)
    blast = "blast"


class BlastTaskChoices(str, Enum):
    megablast = "megablast"
    dc_megablast = "dc-megablast"
    blastn = "blastn"
    blastn_short = "blastn-short"


@app.command()
def run(
    input: str = typer.Option(..., "--input", help="Path to the input FASTA file."),
    output_dir: str = typer.Option(
        DEFAULT_OUTPUT_DIR, "--output-dir", help="Directory to write output files."
    ),
    output_file: str = typer.Option(
        DEFAULT_OUTPUT_FILE,
        "--output-file",
        help="File to write final results. Valid extensions: .csv, .tsv or .json",
    ),
    datasets_dir: str = typer.Option(
        DEFAULT_DATASETS_LOCAL_PATH,
        "--datasets-dir",
        help="Path to local directory containing nextclade datasets.",
    ),
    nextclade_sort_min_score: float = typer.Option(
        DEFAULT_NEXTCLADE_SORT_MIN_SCORE,
        "--ns-min-score",
        help="Nextclade sort min score.",
    ),
    nextclade_sort_min_hits: int = typer.Option(
        DEFAULT_NEXTCLADE_SORT_MIN_HITS,
        "--ns-min-hits",
        help="Nextclade sort min hits.",
    ),
    blast_database: str = typer.Option(
        DEFAULT_BLAST_DATABASE,
        "--blast-database",
        help="Path to local blast database.",
    ),
    blast_database_metadata: str = typer.Option(
        DEFAULT_BLAST_DATABASE_METADATA,
        "--blast-database-metadata",
        help="Path to local blast database metadata.",
    ),
    blast_pident: int = typer.Option(
        DEFAULT_BLAST_IDENTITY_THRESHOLD,
        "--blast-pident",
        help="Identity threshold for BLAST analysis.",
    ),
    blast_evalue: float = typer.Option(
        DEFAULT_BLAST_EVALUE,
        "--blast-evalue",
        help="E-value threshold for BLAST analysis.",
    ),
    blast_qcov: int = typer.Option(
        DEFAULT_BLAST_QCOV,
        "--blast-qcov",
        help="Minimum query coverage per HSP for BLAST analysis.",
    ),
    blast_task: BlastTaskChoices = typer.Option(
        BlastTaskChoices.megablast,
        "--blast-task",
        help="BLAST task type (megablast, dc-megablast, blastn, blastn-short).",
    ),
    config_file_path: Optional[str] = DEFAULT_CONFIG_FILE,
    snk_file_path: Optional[str] = DEFAULT_SNK_FILE,
    cores: int = DEFAULT_CORES,
    verbose: bool = typer.Option(
        DEFAULT_VERBOSE, "--verbose", "-v", help="Show snakemake logs."
    ),
):
    """Identify virus sequences and run nextclade for each virus."""
    snakemake_response = run_analysis.run(
        snk_file=snk_file_path,
        config_file=config_file_path,
        cores=cores,
        sequences_fasta=input,
        output_dir=output_dir,
        output_file=output_file,
        datasets_local_path=datasets_dir,
        nextclade_sort_min_score=nextclade_sort_min_score,
        nextclade_sort_min_hits=nextclade_sort_min_hits,
        blast_database=blast_database,
        blast_database_metadata=blast_database_metadata,
        blast_identity_threshold=blast_pident,
        blast_evalue=blast_evalue,
        blast_qcov=blast_qcov,
        blast_task=blast_task.value,
        verbose=verbose,
    )
    if snakemake_response.status == 200:
        log_multiline(snakemake_response.format_log())
        logger.info("Nextclade run with success.")
    else:
        log_multiline(snakemake_response.format_log())
        logger.error("Failed to run nextclade.")
