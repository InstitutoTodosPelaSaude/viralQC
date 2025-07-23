from snakemake import snakemake
from viralqa.core.utils import capture_snakemake_log
from viralqa.core.errors import SnakemakeExecutionFailed


class GetNextcladeDatasets:
    def __init__(self):
        pass

    @capture_snakemake_log
    def get_public_dataset(
        self,
        snk_file: str,
        config_file: str,
        cores: int,
    ) -> str:
        successful = snakemake(
            snk_file, configfiles=[config_file], cores=cores, targets=["all"]
        )
        if not successful:
            raise SnakemakeExecutionFailed(snk_file, "")
        return "Nextclade public datasets recovered."
