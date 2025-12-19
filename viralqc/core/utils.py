import re
import subprocess
import tempfile
import yaml
import os
from typing import Tuple, Optional
from pathlib import Path
from viralqc.core.models import SnakemakeResponse, RunStatus


def _get_log_path_from_workdir(workdir: str) -> Tuple[str, Optional[str]]:
    """
    Find the most recent Snakemake log file in the workdir.
    Returns (log_path, run_id) tuple.
    """
    if not workdir:
        return "This execution has no working directory.", None

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
    Run Snakemake via subprocess instead of using the Python API.
    This ensures better isolation between concurrent runs.

    Keyword arguments:
        snk_file -- .snk snakemake file path
        config_file -- .yaml config file path
        cores -- number of cores used to run snakemake
        config -- dictionary of config parameters
        workdir -- working directory for snakemake
        verbose -- whether to show verbose output
    """
    cmd = ["snakemake", "-s", snk_file, "-c", str(cores), "all"]
    merged_config = {}

    if config_file:
        config_file_list = (
            config_file if isinstance(config_file, list) else [config_file]
        )
        for cf in config_file_list:
            with open(str(cf), "r") as f:
                user_config = yaml.safe_load(f) or {}
                merged_config.update(user_config)

    if config:
        for key, value in config.items():
            if value is None:
                continue
            if hasattr(value, "__fspath__"):
                merged_config[key] = str(value)
            else:
                merged_config[key] = value

    temp_config_file = None
    if merged_config:
        temp_config_file = tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False
        )
        yaml.dump(merged_config, temp_config_file)
        temp_config_file.flush()
        temp_config_file.close()

        cmd.extend(["--configfile", temp_config_file.name])

    if workdir:
        cmd.extend(["--directory", workdir])

    if not verbose:
        cmd.append("--quiet")

    try:
        if verbose:
            result = subprocess.run(cmd, text=True, check=False)
            captured_output = ""
        else:
            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            captured_output = ""
            if result.stdout:
                captured_output += result.stdout
            if result.stderr:
                captured_output += "\n" + result.stderr

        successful = result.returncode == 0

    except Exception as e:
        successful = False
        captured_output = f"Exception during Snakemake execution: {str(e)}"
    finally:
        if temp_config_file and os.path.exists(temp_config_file.name):
            os.unlink(temp_config_file.name)

    log_path, run_id = _get_log_path_from_workdir(workdir)

    results_path = None
    if config:
        output_dir = config.get("output_dir", "")
        output_file = config.get("output_file", "results.json")
        if output_dir and output_file:
            results_path = f"{output_dir}/outputs/{output_file}"

    if successful:
        return SnakemakeResponse(
            run_id=run_id,
            status=RunStatus.SUCCESS,
            log_path=log_path,
            results_path=results_path,
            captured_output=captured_output,
        )
    else:
        return SnakemakeResponse(
            run_id=run_id,
            status=RunStatus.FAIL,
            log_path=log_path,
            results_path=results_path,
            captured_output=captured_output,
        )
