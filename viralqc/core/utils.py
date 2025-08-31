import io, contextlib, re
from typing import Tuple, Optional, List
from snakemake import snakemake
from viralqc.core.errors import SnakemakeExecutionFailed
from viralqc.core.models import SnakemakeResponse, RunStatus


def _get_log_and_run_id_from_log(log_lines: str) -> Tuple[str, Optional[str]]:
    last_line = log_lines.strip().split("\n")[-1]
    match = re.search(r"\d{4}-\d{2}-\d{2}T\d{6}\.\d+", last_line)
    if "Complete log" in last_line:
        log_path = re.sub("Complete log: ", "", last_line)
    else:
        log_path = "This execution has no log file."
    run_id = match.group() if match else None
    return log_path, run_id


def run_snakemake(
    snk_file: str,
    config_file: Optional[List[str]] = None,
    cores: int = 1,
    config: dict = None,
) -> SnakemakeResponse:
    """
    The snakemake module has runtime logic that must be handled with viralQA
    modularization patterns, including:
        - returns only a Boolean indicating whether the flow ran successfully or not.
        - all logs are output as stderr on the console.

    Therefore, this function handles this.

    Keyword arguments:
        snk_file -- .snk snakemake file path
        config_file -- .yml or .json snakemake config file path
        cores -- number of cores used to run snakemake
    """
    stdout_buf = io.StringIO()
    with contextlib.redirect_stderr(stdout_buf):
        successful = snakemake(
            snk_file,
            config=config,
            configfiles=config_file,
            cores=cores,
            targets=["all"],
        )
        stdout = stdout_buf.getvalue()
        log_path, run_id = _get_log_and_run_id_from_log(stdout)

        try:
            if successful:
                return SnakemakeResponse(
                    run_id=run_id,
                    script_path=snk_file,
                    config_path=config_file,
                    status=RunStatus.SUCCESS.value,
                    log_path=log_path,
                    log_content=stdout,
                )
            else:
                return SnakemakeResponse(
                    run_id=run_id,
                    script_path=snk_file,
                    config_path=config_file,
                    status=RunStatus.FAIL.value,
                    log_path=log_path,
                    log_content=stdout,
                    error=SnakemakeExecutionFailed(snk_file),
                )
        except:
            return SnakemakeResponse(
                run_id=run_id,
                script_path=snk_file,
                config_path=config_file,
                status=RunStatus.FAIL.value,
                log_path=log_path,
                log_content=stdout,
                error=Exception,
            )
