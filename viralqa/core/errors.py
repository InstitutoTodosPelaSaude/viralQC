class SnakemakeExecutionFailed(Exception):
    """
    Raised when a snakemake() invocation fails.
    """

    def __init__(self, script):
        msg = f"Snakemake run of {script} failed"
        super().__init__(msg)
