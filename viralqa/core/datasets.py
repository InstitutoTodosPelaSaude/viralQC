from viralqa.core.utils import run_snakemake
from viralqa.core.models import SnakemakeResponse


class GetNextcladeDatasets:
    def __init__(self):
        pass

    def get_public_dataset(
        self,
        datasets_dir: str,
        snk_file: str,
        config_file: str,
        cores: int,
    ) -> SnakemakeResponse:
        config = {"datasets_dir": datasets_dir}
        snakemake_response = run_snakemake(snk_file, [config_file], cores, config)
        return snakemake_response


class GetBlastDatabase:
    def __init__(self):
        pass

    def get_database(
        self,
        datasets: str,
        output_dir: str,
        snk_file: str,
        cores: int,
    ) -> SnakemakeResponse:
        config = {"datasets": datasets, "output_dir": output_dir}

        snakemake_response = run_snakemake(snk_file, None, cores, config)
        return snakemake_response
