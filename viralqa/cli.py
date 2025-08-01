import typer, logging, colorlog
from typing import Optional
from viralqa.core.datasets import GetNextcladeDatasets
from viralqa import GET_NC_PUBLIC_DATASETS_SNK_PATH, GET_NC_PUBLIC_DATASETS_CONFIG_PATH

# core config
get_nc_datasets = GetNextcladeDatasets()

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
logger.setLevel(logging.DEBUG)
logger.addHandler(handler)

# cli config
app = typer.Typer()


@app.command()
def get_nextclade_datasets(
    snk_file_path: Optional[str] = GET_NC_PUBLIC_DATASETS_SNK_PATH,
    config_file_path: Optional[str] = GET_NC_PUBLIC_DATASETS_CONFIG_PATH,
    cores: int = 1,
):
    """Get Nextclade virus dataset"""
    snakemake_response = get_nc_datasets.get_public_dataset(
        snk_file=snk_file_path, config_file=config_file_path, cores=cores
    )
    if snakemake_response.status == 200:
        logger.info(snakemake_response.format_log())
        logger.info("Nextclade public datasets successfully retrieved.")
    else:
        logger.error(snakemake_response.format_log())
        logger.error("Failed to retrieve Nextclade public datasets.")


@app.command()
def get_custom_datasets(cores: int = 1):
    """Get custom virus dataset"""
    print("In progress")


if __name__ == "__main__":
    app()
