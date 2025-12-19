from viralqc.core.utils import run_snakemake
from viralqc.core.models import SnakemakeResponse
from viralqc.core.defaults import DEFAULT_VERBOSE


class GetNextcladeDatasets:
    def __init__(self):
        pass

    def get_public_dataset(
        self,
        datasets_dir: str,
        snk_file: str,
        config_file: str,
        cores: int,
        verbose: bool = DEFAULT_VERBOSE,
    ) -> SnakemakeResponse:
        config = {"datasets_dir": datasets_dir}
        snakemake_response = run_snakemake(
            snk_file=snk_file,
            config_file=[config_file],
            cores=cores,
            config=config,
            workdir=None,
            verbose=verbose,
        )
        return snakemake_response


class GetBlastDatabase:
    def __init__(self):
        pass

    def get_database(
        self,
        output_dir: str,
        release_date: str | None,
        snk_file: str,
        cores: int,
        verbose: bool = DEFAULT_VERBOSE,
    ) -> SnakemakeResponse:
        config = {
            "output_dir": output_dir,
            "release_date": release_date,
        }

        snakemake_response = run_snakemake(
            snk_file=snk_file,
            config_file=None,
            cores=cores,
            config=config,
            workdir=None,
            verbose=verbose,
        )
        return snakemake_response
