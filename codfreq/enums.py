"""Enumeration types used across CodFreq CLI applications."""

from enum import Enum


class LogFormat(str, Enum):
    """Supported output formats for command-line logging."""

    text = "text"
    json = "json"


class Program(str, Enum):
    """Available alignment programs."""

    minimap2 = "minimap2"
    bowtie2 = "bowtie2"
