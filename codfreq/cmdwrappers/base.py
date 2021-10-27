import click  # type: ignore
from subprocess import Popen, PIPE

REFINIT_FUNCTIONS = {}
ALIGN_FUNCTIONS = {}
AUTOREMOVE_CONTAINERS = False


def execute(command):
    proc = Popen(command, stdout=PIPE, stderr=PIPE, encoding='U8')
    out, err = proc.communicate()
    raise_on_proc_error(proc, err)
    return out, err


def raise_on_proc_error(proc, err):
    if proc.returncode:
        click.echo(err, err=True)
        raise click.Abort()


def refinit_func(name):

    def wrapper(func):
        REFINIT_FUNCTIONS[name] = func
        return func

    return wrapper


def align_func(name):

    def wrapper(func):
        ALIGN_FUNCTIONS[name] = func
        return func

    return wrapper


def get_programs():
    return sorted(ALIGN_FUNCTIONS.keys())


def get_refinit(name):
    return REFINIT_FUNCTIONS[name]


def get_align(name):
    return ALIGN_FUNCTIONS[name]
