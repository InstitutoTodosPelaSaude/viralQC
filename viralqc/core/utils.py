import re
from typing import Tuple, Optional
from pathlib import Path
from snakemake import snakemake
from viralqc.core.models import SnakemakeResponse, RunStatus


def _get_log_path_from_workdir(workdir: str) -> Tuple[str, Optional[str]]:
    """
    Find the most recent Snakemake log file in the workdir.
    Returns (log_path, run_id) tuple.
    """
    if not workdir:
        return "This execution has no log file.", None

    log_dir = Path(workdir) / ".snakemake" / "log"

    if not log_dir.exists():
        return "This execution has no log file.", None

    log_files = list(log_dir.glob("*.snakemake.log"))

    if not log_files:
        return "This execution has no log file.", None
    most_recent_log = max(log_files, key=lambda p: p.stat().st_mtime)
    log_path = str(most_recent_log)

    match = re.search(r"(\d{4}-\d{2}-\d{2}T\d{6}\.\d+)", most_recent_log.name)
    run_id = match.group(1) if match else None

    return log_path, run_id


def run_snakemake(
    snk_file: str,
    config_file: Path | None = None,
    cores: int = 1,
    config: dict = None,
    workdir: str = None,
    verbose: bool = False,
) -> SnakemakeResponse:
    """
    The snakemake module has runtime logic that must be handled with viralQA
    modularization patterns, including:
        - returns only a Boolean indicating whether the flow ran successfully or not.
        - all logs are output as stderr on the console.

    Therefore, this function handles this.

    Keyword arguments:
        snk_file -- .snk snakemake file path
        config_file -- .yaml config file path
        cores -- number of cores used to run snakemake
    """
    successful = snakemake(
        snk_file,
        config=config,
        configfiles=config_file,
        cores=cores,
        targets=["all"],
        workdir=workdir,
        quiet=not verbose,
    )

    log_path, run_id = _get_log_path_from_workdir(workdir)

    results_path = None
    if config:
        output_dir = config.get("output_dir", "")
        output_file = config.get("output_file", "results.json")
        if output_dir and output_file:
            results_path = f"{output_dir}/outputs/{output_file}"

    try:
        if successful:
            return SnakemakeResponse(
                run_id=run_id,
                status=RunStatus.SUCCESS,
                log_path=log_path,
                results_path=results_path,
                captured_output="",
            )
        else:
            return SnakemakeResponse(
                run_id=run_id,
                status=RunStatus.FAIL,
                log_path=log_path,
                results_path=results_path,
                captured_output="",
            )
    except Exception as e:
        return SnakemakeResponse(
            run_id=run_id,
            status=RunStatus.FAIL,
            log_path=log_path,
            results_path=results_path,
            captured_output=f"Exception during Snakemake execution: {str(e)}",
        )
