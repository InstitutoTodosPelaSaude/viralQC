class SnakemakeExecutionFailed(Exception):
    """
    Raised when a snakemake() invocation fails.
    Carries both the snakemake script name and the full stdout+stderr.
    """

    def __init__(self, script: str, log: str):
        msg = f"Snakemake run of {script} failed.\n\n--- log start ---\n{log}\n--- log end ---"
        super().__init__(msg)
        self.script = script
        self.log = log
