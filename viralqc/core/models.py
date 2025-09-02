from dataclasses import dataclass, field
from enum import IntEnum
from typing import Optional, List, Union
from viralqc.core.errors import SnakemakeExecutionFailed


class RunStatus(IntEnum):
    SUCCESS = 200
    INPROGRESS = 202
    FAIL = 500


@dataclass
class SnakemakeResponse:
    """
    Represents the result of a Snakemake workflow execution.

    Attributes:
        run_id (str): Unique identifier for the execution run, derived from a timestamp in snakemake log.
        script_path (str): Path to the Snakemake workflow script (.smk) used for the run.
        config_path (str): Path to the configuration file (.yml, .json) used during the execution.
        status (RunStatus): Execution status indicating success, in progress or failure.
        log (str): Combined stdout and stderr output captured during the run.
        error (Optional[List[Union[SnakemakeExecutionFailed, Exception]]]):
            A list of errors raised during execution, including custom or generic exceptions.
            Defaults to None if the execution was successful.
    """

    run_id: str
    script_path: str
    config_path: str
    status: RunStatus
    log_path: str
    log_content: str
    error: Optional[List[Union[SnakemakeExecutionFailed, Exception]]] = field(
        default=None
    )

    def format_log(self) -> str:
        return self.log_content.replace("\\n", "\n")
