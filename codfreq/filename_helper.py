import os
import re
from typing import Optional, Tuple

from .codfreq_types import FASTQFileName, Header

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


def name_bamfile(
    name: str,
    ref_name: Header,
    is_trimmed: bool = True
) -> str:
    return (
        '{}.{}.bam' if is_trimmed else '{}.{}.orig.bam'
    ).format(name, ref_name)


def name_codfreq(name: str) -> str:
    return '{}.codfreq'.format(name)


def name_segfreq(name: str, refname: Header) -> str:
    return '{}.{}.segfreq'.format(name, refname)


def name_nucfreq(name: str) -> str:
    return '{}.nucfreq'.format(name)


def name_binucfreq(name: str) -> str:
    return '{}.binucfreq'.format(name)


def name_patterns(name: str, frag_name: str) -> str:
    return '{}.{}-patterns.fasta'.format(name, frag_name)


def name_consensus(dirname: str, gene: str, level: float) -> str:
    return os.path.join(
        dirname,
        'consensus-{}-{:g}.fasta'.format(gene, level * 100)
    )


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
