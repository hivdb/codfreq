#! /usr/bin/env python

import pysam  # type: ignore
from pysam import AlignedSegment  # type: ignore
from collections import OrderedDict

from .codfreq_types import Header

__all__ = ['PairedReads', 'iter_paired_reads']


PairedReads = tuple[str, list[pysam.AlignedSegment]]


def iter_paired_reads(
    samfile: str
) -> list[PairedReads]:
    idx: int
    name: Header | None
    read: AlignedSegment
    paired_reads: dict[str, list[AlignedSegment]] = OrderedDict()
    with pysam.AlignmentFile(samfile, 'rb') as fp:
        for idx, read in enumerate(fp.fetch()):
            name = read.query_name
            if name is None:
                continue
            paired_reads.setdefault(name, []).append(read)
    return list(paired_reads.items())
