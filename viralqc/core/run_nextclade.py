from viralqc.core.utils import run_snakemake
from viralqc.core.models import SnakemakeResponse
from viralqc.core.errors import InvalidOutputFormat


class RunNextclade:
    def __init__(self):
        pass

    def _get_output_format(self, output_file: str) -> str | InvalidOutputFormat:
        file_extension = output_file.split(".")[-1]
        if file_extension not in ["csv", "tsv", "json"]:
            raise InvalidOutputFormat(file_extension)
        return file_extension

    def run(
        self,
        snk_file: str,
        config_file: str,
        cores: int,
        sequences_fasta: str,
        output_dir: str,
        output_file: str,
        datasets_local_path: str,
        nextclade_sort_min_score: float,
        nextclade_sort_min_hits: int,
        blast_database: str,
        blast_identity_threshold: float,
    ) -> SnakemakeResponse:
        output_format = self._get_output_format(output_file)
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
            "blast_identity_threshold": blast_identity_threshold,
        }

        snakemake_response = run_snakemake(snk_file, [config_file], cores, config)
        return snakemake_response
