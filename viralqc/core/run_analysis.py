from viralqc.core.utils import run_snakemake
from viralqc.core.models import SnakemakeResponse
from viralqc.core.errors import InvalidOutputFormat
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
from pathlib import Path


class RunAnalysis:
    def __init__(self):
        pass

    def _get_output_format(self, output_file: str) -> str | InvalidOutputFormat:
        file_extension = output_file.split(".")[-1]
        if file_extension not in ["csv", "tsv", "json"]:
            raise InvalidOutputFormat(file_extension)
        return file_extension

    def run(
        self,
        sequences_fasta: Path,
        snk_file: str = DEFAULT_SNK_FILE,
        config_file: str = DEFAULT_CONFIG_FILE,
        cores: int = DEFAULT_CORES,
        output_dir: str = DEFAULT_OUTPUT_DIR,
        output_file: Path = DEFAULT_OUTPUT_FILE,
        datasets_local_path: Path = DEFAULT_DATASETS_LOCAL_PATH,
        nextclade_sort_min_score: float = DEFAULT_NEXTCLADE_SORT_MIN_SCORE,
        nextclade_sort_min_hits: int = DEFAULT_NEXTCLADE_SORT_MIN_HITS,
        blast_database: Path = DEFAULT_BLAST_DATABASE,
        blast_database_metadata: Path = DEFAULT_BLAST_DATABASE_METADATA,
        blast_identity_threshold: float = DEFAULT_BLAST_IDENTITY_THRESHOLD,
        blast_evalue: float = DEFAULT_BLAST_EVALUE,
        blast_qcov: float = DEFAULT_BLAST_QCOV,
        blast_task: str = DEFAULT_BLAST_TASK,
        verbose: bool = DEFAULT_VERBOSE,
    ) -> SnakemakeResponse:
        output_format = self._get_output_format(str(output_file))
        config = {
            "sequences_fasta": sequences_fasta,
            "output_dir": output_dir,
            "output_file": output_file,
            "output_format": output_format,
            "config_file": config_file,
            "datasets_local_path": datasets_local_path,
            "threads": cores,
            "nextclade_sort_min_score": nextclade_sort_min_score,
            "nextclade_sort_min_hits": nextclade_sort_min_hits,
            "blast_database": blast_database,
            "blast_database_metadata": blast_database_metadata,
            "blast_identity_threshold": blast_identity_threshold,
            "blast_evalue": blast_evalue,
            "blast_qcov": blast_qcov,
            "blast_task": blast_task,
        }

        snakemake_response = run_snakemake(
            snk_file, [config_file], cores, config, verbose
        )
        return snakemake_response
