from fastapi import FastAPI, Query, HTTPException
from viralqc.core.datasets import GetNextcladeDatasets
from viralqc import GET_NC_PUBLIC_DATASETS_SNK_PATH, DATASETS_CONFIG_PATH

app = FastAPI(
    title="ViralQC Example API",
    description="A demo REST API for the viralQA.",
)

get_nc_datasets = GetNextcladeDatasets()


@app.get("/")
def root():
    return {"message": "Welcome to ViralQC API!"}


@app.get("/get_nextclade_datasets")
def get_nextclade_datasets(cores: int = Query(...)):
    """Get nextclade datasets"""
    snakemake_response = get_nc_datasets.get_public_dataset(
        datasets_dir="datasets",
        snk_file=GET_NC_PUBLIC_DATASETS_SNK_PATH,
        config_file=DATASETS_CONFIG_PATH,
        cores=cores,
    )
    if snakemake_response.status == 200:
        return {"result": snakemake_response.format_log()}
    else:
        raise HTTPException(
            status_code=500,
            detail=f"Snakemake execution error: {str(snakemake_response.format_log())}",
        )


def start():
    import uvicorn

    uvicorn.run("viralqc.api:app", host="127.0.0.1", port=8000, reload=True)
