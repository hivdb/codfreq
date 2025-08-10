"""Enumeration types used across CodFreq CLI applications."""

from enum import Enum

from .cmdwrappers import get_programs


class LogFormat(str, Enum):
    """Supported output formats for command-line logging."""

    text = "text"
    json = "json"


Program = Enum(  # type: ignore[misc]
    "Program",
    {name.upper(): name for name in get_programs()},
    type=str,
    module=__name__,
)
Program.__doc__ = "Available alignment programs."

