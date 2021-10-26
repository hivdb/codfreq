import os
import re
from typing import Optional, Tuple

from .codfreq_types import FASTQFileName, GeneText

FnPair = Tuple[FASTQFileName, Optional[FASTQFileName]]
Pattern = Tuple[
    str,  # delimiter
    int,  # diffoffset
    int,  # pos_paired_marker
    int,  # reverse
]


def suggest_pair_name(
    fnpair: FnPair,
    pattern: Pattern
) -> str:
    filename, _ = fnpair
    delimiter, offset, _, reverse = pattern
    dirpath, filename = os.path.split(filename)
    samfile = re.split(r'(?i)\.fastq(?:.gz)?$', filename)[0]
    if reverse == -1:
        return os.path.join(
            dirpath, samfile)
    samfile = samfile.split(delimiter)
    if reverse:
        samfile.reverse()
    samfile = samfile[:offset] + samfile[offset + 1:]
    if reverse:
        samfile.reverse()
    return os.path.join(
        dirpath, delimiter.join(samfile))


def name_file(
    fnpair: FnPair,
    pattern: Pattern,
    suffix: str
) -> str:
    return suggest_pair_name(fnpair, pattern) + suffix


def name_bamfile(name: str, gene: GeneText) -> str:
    return '{}.{}.bam'.format(name, gene)


def name_codfreq(name: str) -> str:
    return '{}.codfreq'.format(name)


def replace_ext(
        filename: str,
        toext: str,
        fromext: Optional[str] = None,
        name_only: bool = False
) -> str:
    if name_only:
        filename = os.path.split(filename)[-1]
    if fromext:
        return filename[-len(fromext):] + toext
    else:
        return os.path.splitext(filename)[0] + toext
