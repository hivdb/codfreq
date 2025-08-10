import sys
from subprocess import Popen, PIPE
from collections.abc import Callable

import rich
import typer

REFINIT_FUNCTIONS = {}
ALIGN_FUNCTIONS = {}
AUTOREMOVE_CONTAINERS = False


def execute(command: list[str]) -> tuple[str, str]:
    """Execute a subprocess and capture its output.

    :param command: Command and arguments to run.
    :returns: Tuple of standard output and error text.
    :raises typer.Abort: If the command exits with a non-zero status.
    """
    proc = Popen(command, stdout=PIPE, stderr=PIPE, encoding='U8')
    out, err = proc.communicate()
    raise_on_proc_error(proc, err)
    return out, err


def raise_on_proc_error(proc: Popen, err: str) -> None:
    """Abort the program if the subprocess returned an error.

    :param proc: Completed subprocess instance.
    :param err: Captured standard error output.
    :returns: None
    :raises typer.Abort: If the subprocess had a non-zero return code.
    """
    if proc.returncode:
        rich.print(err, file=sys.stderr)
        raise typer.Abort()


def refinit_func(name: str) -> Callable:

    def wrapper(func: Callable) -> Callable:
        REFINIT_FUNCTIONS[name] = func
        return func

    return wrapper


def align_func(name: str) -> Callable:

    def wrapper(func: Callable) -> Callable:
        ALIGN_FUNCTIONS[name] = func
        return func

    return wrapper


def get_programs() -> list[str]:
    return sorted(ALIGN_FUNCTIONS.keys())


def get_refinit(name: str) -> Callable:
    return REFINIT_FUNCTIONS[name]


def get_align(name: str) -> Callable:
    return ALIGN_FUNCTIONS[name]
