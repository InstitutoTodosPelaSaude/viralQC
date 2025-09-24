from fastapi import FastAPI, Query, HTTPException, UploadFile, File
from viralqc.core.datasets import GetNextcladeDatasets, GetBlastDatabase
from pathlib import Path
from uuid import uuid4
from json import load
from viralqc.core.run_nextclade import RunNextclade
from viralqc import (
    DATASETS_CONFIG_PATH,
    GET_NC_PUBLIC_DATASETS_SNK_PATH,
    GET_BLAST_DB_SNK_PATH,
    RUN_NEXTCLADE_SNK_PATH,
)

app = FastAPI(
    title="ViralQC Example API",
    description="A REST API for the viralQC.",
)

get_nc_datasets = GetNextcladeDatasets()
get_blast_db = GetBlastDatabase()
run_nextclade = RunNextclade()


def _get_tmp_dir_uuid() -> Path:
    tmp_dir = Path("/tmp/vqc") / str(uuid4())
    tmp_dir.mkdir(parents=True, exist_ok=True)
    return tmp_dir


def _save_input_file(content: str, dir: Path) -> Path:
    file_path = dir / "input.fasta"
    with open(file_path, "wb") as f:
        f.write(content)
    return file_path


@app.get("/")
def root():
    return {"message": "Welcome to ViralQC API!"}


@app.post("/get_nextclade_datasets")
def get_nextclade_datasets(cores: int = Query(1, description="Number of cores to use")):
    """Create local nextclade datasets"""
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


@app.post("/get_blast_database")
def get_blast_database(cores: int = Query(1, description="Number of cores to use")):
    """Create BLAST database based on ncbi viruses refseq genomes"""
    snakemake_response = get_blast_db.get_database(
        output_dir="datasets",
        snk_file=GET_BLAST_DB_SNK_PATH,
        cores=cores,
    )
    if snakemake_response.status == 200:
        return {"result": snakemake_response.format_log()}
    else:
        raise HTTPException(
            status_code=500,
            detail=f"Snakemake execution error: {str(snakemake_response.format_log())}",
        )


@app.post("/run")
async def run(
    sequences_fasta: UploadFile = File(...),
    cores: int = Query(1, description="Number of cores to use"),
    nextclade_sort_min_score: float = Query(
        0.1, description="Nextclade sort min score."
    ),
    nextclade_sort_min_hits: int = Query(10, description="Nextclade sort min hits."),
    identity_threshold: float = Query(
        0.9, description="Minimum percentage threshold to BLAST consider as hit."
    ),
):
    """Run viralQC given an input fasta file."""
    output_directory = _get_tmp_dir_uuid()
    file_binary_content = await sequences_fasta.read()
    input_file = _save_input_file(file_binary_content, output_directory)
    snakemake_response = run_nextclade.run(
        snk_file=RUN_NEXTCLADE_SNK_PATH,
        config_file=DATASETS_CONFIG_PATH,
        cores=cores,
        sequences_fasta=input_file,
        output_dir=output_directory,
        output_file="results.json",
        datasets_local_path="datasets",
        nextclade_sort_min_score=nextclade_sort_min_score,
        nextclade_sort_min_hits=nextclade_sort_min_hits,
        blast_database="datasets/blast.fasta",
        blast_identity_threshold=identity_threshold,
    )
    with open(f"{output_directory}/results.json", "r") as f:
        results_data = load(f)
    if snakemake_response.status == 200:
        return {
            "log": snakemake_response.format_log(),
            "tmp_results_path": output_directory,
            "sequences_to_submmit": f"{output_directory}/sequences_target_regions.fasta",
            "results": results_data,
        }
    else:
        raise HTTPException(
            status_code=500,
            detail=f"Snakemake execution error: {str(snakemake_response.format_log())}",
        )


def start():
    import uvicorn

    uvicorn.run("viralqc.api:app", host="127.0.0.1", port=8000, reload=True)
