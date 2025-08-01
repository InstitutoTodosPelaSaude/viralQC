from viralqa.core.utils import run_snakemake
from viralqa.core.models import SnakemakeResponse


class GetNextcladeDatasets:
    def __init__(self):
        pass

    def get_public_dataset(
        self,
        snk_file: str,
        config_file: str,
        cores: int,
    ) -> SnakemakeResponse:
        snakemake_response = run_snakemake(snk_file, config_file, cores)
        return snakemake_response
