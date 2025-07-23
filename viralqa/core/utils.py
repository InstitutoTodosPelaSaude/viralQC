import functools
import io
import contextlib
from viralqa.core.errors import SnakemakeExecutionFailed


def capture_snakemake_log(func):
    """
    Captures stdout+stderr during the wrapped function’s execution.
    If SnakemakeExecutionFailed is raised inside, re-raise it with the captured log.
    Also returns or prints the log when successful.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        stdout_buf = io.StringIO()
        stderr_buf = io.StringIO()

        with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(
            stderr_buf
        ):
            try:
                result = func(*args, **kwargs)
            except SnakemakeExecutionFailed as e:
                out = stdout_buf.getvalue()
                err = stderr_buf.getvalue()
                full_log = (out + "\n" + err).strip()
                raise SnakemakeExecutionFailed(e.script, full_log) from e

        # If no error print the snakemake log
        out = stdout_buf.getvalue()
        err = stderr_buf.getvalue()
        print(out + "\n" + err)

        return result

    return wrapper
