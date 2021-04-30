#! /usr/bin/env python

import pysam
from collections import OrderedDict


def iter_paired_reads(samfile):
    paired_reads = OrderedDict()
    with pysam.AlignmentFile(samfile, 'rb') as samfile:
        for idx, read in enumerate(samfile.fetch()):
            name = read.query_name
            paired_reads.setdefault(name, []).append(read)
    return paired_reads.items()
