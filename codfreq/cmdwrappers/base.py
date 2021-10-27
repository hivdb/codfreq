import click  # type: ignore
from subprocess import Popen, PIPE
from typing import List, Tuple, Callable

REFINIT_FUNCTIONS = {}
ALIGN_FUNCTIONS = {}
AUTOREMOVE_CONTAINERS = False


def execute(command: List[str]) -> Tuple[str, str]:
    proc = Popen(command, stdout=PIPE, stderr=PIPE, encoding='U8')
    out, err = proc.communicate()
    raise_on_proc_error(proc, err)
    return out, err


def raise_on_proc_error(proc: Popen, err: str) -> None:
    if proc.returncode:
        click.echo(err, err=True)
        raise click.Abort()


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


def get_programs() -> List[str]:
    return sorted(ALIGN_FUNCTIONS.keys())


def get_refinit(name: str) -> Callable:
    return REFINIT_FUNCTIONS[name]


def get_align(name: str) -> Callable:
    return ALIGN_FUNCTIONS[name]
