import io, contextlib, re, sys
from typing import Tuple, Optional
from pathlib import Path
from snakemake import snakemake
from viralqc.core.models import SnakemakeResponse, RunStatus


class Tee(object):
    """
    Redirects output to multiple files, the flush method is called for each file
    in order to ensure that the output is written to all files and that the
    output is not buffered, so in the context of this code the content is directly
    written to the log files and printed to the console (if verbose=True).
    """

    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()

    def flush(self):
        for f in self.files:
            f.flush()


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
    config_file: Path | None = None,
    cores: int = 1,
    config: dict = None,
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
    stdout_buf = io.StringIO()

    if verbose:
        ctx = contextlib.redirect_stderr(Tee(sys.stderr, stdout_buf))
    else:
        ctx = contextlib.redirect_stderr(stdout_buf)

    with ctx:
        successful = snakemake(
            snk_file,
            config=config,
            configfiles=config_file,
            cores=cores,
            targets=["all"],
        )
        stdout = stdout_buf.getvalue()
        log_path, run_id = _get_log_and_run_id_from_log(stdout)

        # Construct results_path from config
        results_path = None
        if config:
            output_dir = config.get("output_dir", "")
            output_file = config.get("output_file", "results.json")
            if output_dir and output_file:
                results_path = f"{output_dir}/{output_file}"

        try:
            if successful:
                return SnakemakeResponse(
                    run_id=run_id,
                    status=RunStatus.SUCCESS,
                    log_path=log_path,
                    results_path=results_path,
                    captured_output=stdout,
                )
            else:
                return SnakemakeResponse(
                    run_id=run_id,
                    status=RunStatus.FAIL,
                    log_path=log_path,
                    results_path=results_path,
                    captured_output=stdout,
                )
        except Exception as e:
            return SnakemakeResponse(
                run_id=run_id,
                status=RunStatus.FAIL,
                log_path=log_path,
                results_path=results_path,
                captured_output=stdout_buf.getvalue(),
            )
