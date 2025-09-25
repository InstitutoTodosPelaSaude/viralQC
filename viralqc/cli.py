import typer, logging, colorlog
from typing import Optional
from enum import Enum
from viralqc.core.datasets import GetNextcladeDatasets, GetBlastDatabase
from viralqc.core.run_nextclade import RunNextclade
from viralqc import (
    DATASETS_CONFIG_PATH,
    GET_NC_PUBLIC_DATASETS_SNK_PATH,
    GET_BLAST_DB_SNK_PATH,
    RUN_NEXTCLADE_SNK_PATH,
)

# core config
get_nc_datasets = GetNextcladeDatasets()
get_blast_db = GetBlastDatabase()
run_nextclade = RunNextclade()

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
    config_file_path: Optional[str] = DATASETS_CONFIG_PATH,
    cores: int = 1,
):
    """Get Nextclade virus datasets"""
    snakemake_response = get_nc_datasets.get_public_dataset(
        datasets_dir=datasets_dir,
        snk_file=snk_file_path,
        config_file=config_file_path,
        cores=cores,
    )
    if snakemake_response.status == 200:
        log_multiline(snakemake_response.format_log())
        logger.info("Nextclade public datasets successfully retrieved.")
    else:
        logger.error(snakemake_response.format_log())
        logger.error("Failed to retrieve Nextclade public datasets.")


@app.command()
def get_custom_datasets(cores: int = 1):
    """Get custom virus datasets"""
    print("In progress")


@app.command()
def get_blast_database(
    output_dir: str = typer.Option(
        "datasets",
        "--output-dir",
        help="Path to store BLAST database.",
    ),
    snk_file_path: Optional[str] = GET_BLAST_DB_SNK_PATH,
    cores: int = 1,
):
    """Create BLAST database based on ncbi viruses refseq genomes"""
    snakemake_response = get_blast_db.get_database(
        output_dir=output_dir,
        snk_file=snk_file_path,
        cores=cores,
    )
    if snakemake_response.status == 200:
        log_multiline(snakemake_response.format_log())
        logger.info("BLAST database created.")
    else:
        logger.error(snakemake_response.format_log())
        logger.error("Failed to create BLAST database.")


class SortChoices(str, Enum):
    nextclade = ("nextclade",)
    blast = "blast"


@app.command()
def run_from_fasta(
    sequences_fasta: str = typer.Option(
        ..., "--sequences-fasta", help="Path to the input FASTA file."
    ),
    output_dir: str = typer.Option(
        "output", "--output-dir", help="Directory to write output files."
    ),
    output_file: str = typer.Option(
        "results.tsv",
        "--output-file",
        help="File to write final results. Valid extensions: .csv, .tsv or .json",
    ),
    datasets_dir: str = typer.Option(
        "datasets",
        "--datasets-dir",
        help="Path to local directory containing nextclade datasets.",
    ),
    nextclade_sort_min_score: float = typer.Option(
        0.1,
        "--ns-min-score",
        help="Nextclade sort min score.",
    ),
    nextclade_sort_min_hits: int = typer.Option(
        10,
        "--ns-min-hits",
        help="Nextclade sort min hits.",
    ),
    blast_database: str = typer.Option(
        "datasets/blast.fasta",
        "--blast-database",
        help="Path to local blast database.",
    ),
    identity_threshold: str = typer.Option(
        0.90,
        "--identity-threshold",
        help="Identity threshold for BLAST analysis.",
    ),
    config_file_path: Optional[str] = DATASETS_CONFIG_PATH,
    snk_file_path: Optional[str] = RUN_NEXTCLADE_SNK_PATH,
    cores: int = 1,
):
    """Split sequences by viruses and run nextclade for each virus."""
    snakemake_response = run_nextclade.run(
        snk_file=snk_file_path,
        config_file=config_file_path,
        cores=cores,
        sequences_fasta=sequences_fasta,
        output_dir=output_dir,
        output_file=output_file,
        datasets_local_path=datasets_dir,
        nextclade_sort_min_score=nextclade_sort_min_score,
        nextclade_sort_min_hits=nextclade_sort_min_hits,
        blast_database=blast_database,
        blast_identity_threshold=identity_threshold,
    )
    if snakemake_response.status == 200:
        log_multiline(snakemake_response.format_log())
        logger.info("Nextclade run with success.")
    else:
        logger.error(snakemake_response.format_log())
        logger.error("Failed to run nextclade.")
