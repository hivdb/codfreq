#! /usr/bin/env python

import pysam  # type: ignore
from pysam import AlignedSegment  # type: ignore
from collections import OrderedDict
from typing import Tuple, List, Dict

from .codfreq_types import Header

__all__ = ['PairedReads', 'iter_paired_reads']


PairedReads = Tuple[str, List[pysam.AlignedSegment]]


def iter_paired_reads(
    samfile: str
) -> List[PairedReads]:
    idx: int
    name: Header
    read: AlignedSegment
    paired_reads: Dict[str, List[AlignedSegment]] = OrderedDict()
    with pysam.AlignmentFile(samfile, 'rb') as fp:
        for idx, read in enumerate(fp.fetch()):
            name = read.query_name
            paired_reads.setdefault(name, []).append(read)
    return list(paired_reads.items())
