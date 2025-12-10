from dataclasses import dataclass
from enum import IntEnum
from typing import Optional
import json


class RunStatus(IntEnum):
    SUCCESS = 200
    FAIL = 500


@dataclass
class SnakemakeResponse:
    """
    Represents the result of a Snakemake workflow execution.

    Attributes:
        run_id (str): Unique identifier for the execution run, derived from a timestamp in snakemake log.
        status (RunStatus): Execution status, which can be 200 (success) or 500 (failure).
        log_path (str): Path to the log file.
        results_path (Optional[str]): Path to the results file.
        captured_output (str): Error captured from the execution.
    """

    run_id: str
    status: RunStatus
    log_path: str
    results_path: Optional[str] = None
    captured_output: str = ""

    def format_log(self) -> str:
        """Returns the log file content formatted."""
        try:
            with open(self.log_path, "r") as f:
                return f.read()
        except FileNotFoundError:
            if self.captured_output:
                return f"Log file not found. Captured Output:\n{self.captured_output}"
            return "Log file not found."
        except Exception as e:
            if self.captured_output:
                return f"Error reading log file: {str(e)}\nCaptured Output:\n{self.captured_output}"
            return f"Error reading log file: {str(e)}"

    def get_results(self) -> dict:
        """Returns the results file content in dictionary format."""
        if not self.results_path:
            return {"message": "No results file generated for this execution."}

        try:
            with open(self.results_path, "r") as f:
                return json.load(f).get("data")
        except FileNotFoundError:
            return {"error": "Results file not found."}
        except json.JSONDecodeError as e:
            return {"error": f"Invalid JSON in results file: {str(e)}"}
        except Exception as e:
            return {"error": f"Error reading results file: {str(e)}"}
