class SnakemakeExecutionFailed(Exception):
    """
    Raised when a snakemake() invocation fails.
    """

    def __init__(self, script):
        msg = f"Snakemake run of {script} failed"
        super().__init__(msg)


class InvalidOutputFormat(Exception):
    """
    Raised when an unexpected output format is parsed.
    """

    def __init__(self, format):
        msg = f"Invalid output format extension: .{format}"
        super().__init__(msg)
