import os
import re

from .codfreq_types import FASTQFileName, Header

FnPair = tuple[FASTQFileName, FASTQFileName | None]
Pattern = tuple[
    str,  # delimiter
    int,  # diffoffset
    int,  # pos_paired_marker
    int,  # reverse
]


def suggest_pair_name(
    fnpair: FnPair,
    pattern: Pattern
) -> str:
    """Suggest a SAM/BAM base name from a pair of FASTQ files."""

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
    """Generate a file name for a fragment pair."""

    return suggest_pair_name(fnpair, pattern) + suffix


def name_bamfile(
    name: str,
    ref_name: Header,
    is_trimmed: bool = True
) -> str:
    """Return the BAM file name for a fragment."""

    return (
        '{}.{}.bam' if is_trimmed else '{}.{}.orig.bam'
    ).format(name, ref_name)


def name_codfreq(name: str) -> str:
    """Return the CodFreq file name for a fragment."""

    return f'{name}.codfreq'


def replace_ext(
        filename: str,
        toext: str,
        fromext: str | None = None,
        name_only: bool = False
) -> str:
    """Replace a file extension.

    :param filename: Original file name.
    :param toext: New extension including leading dot.
    :param fromext: Expected original extension; if provided the last
        characters are replaced without removing the existing extension.
    :param name_only: If ``True`` only the base name is processed.
    :returns: File name with the new extension.
    """

    if name_only:
        filename = os.path.split(filename)[-1]
    if fromext:
        return filename[-len(fromext):] + toext
    else:
        return os.path.splitext(filename)[0] + toext
