from viralqa.core.utils import run_snakemake
from viralqa.core.models import SnakemakeResponse


class RunNextclade:
    def __init__(self):
        pass

    def run(
        self,
        snk_file: str,
        config_file: str,
        cores: int,
        sequences_fasta: str,
        sort_mode: str,
        output_dir: str,
        datasets_local_path: str,
        nextclade_sort_min_score=float,
        nextclade_sort_min_hits=int,
        blast_database=str,
    ) -> SnakemakeResponse:
        config = {
            "sequences_fasta": sequences_fasta,
            "sort_mode": sort_mode,
            "output_dir": output_dir,
            "config_file": config_file,
            "datasets_local_path": datasets_local_path,
            "threads": cores,
            "nextclade_sort_min_score": nextclade_sort_min_score,
            "nextclade_sort_min_hits": nextclade_sort_min_hits,
            "blast_database": blast_database,
        }

        snakemake_response = run_snakemake(snk_file, [config_file], cores, config)
        return snakemake_response
