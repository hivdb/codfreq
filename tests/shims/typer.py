"""Minimal Typer stub for testing without external dependency."""
from typing import Any, Callable


class Abort(Exception):
    """Simulate Typer's Abort exception."""
    pass


class Typer:
    """Simple stand-in for Typer application."""

    def command(self, *args: Any, **kwargs: Any) -> Callable:
        def decorator(func: Callable) -> Callable:
            return func
        return decorator

    def __call__(self, *args: Any, **kwargs: Any) -> None:
        raise SystemExit


def echo(message: Any, err: bool | None = False) -> None:
    print(message)


def Argument(*args: Any, **kwargs: Any) -> Any:
    return kwargs.get("default")


def Option(default: Any = None, *args: Any, **kwargs: Any) -> Any:
    return default
