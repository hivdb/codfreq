#! /usr/bin/env python

import pysam
from collections import OrderedDict


def iter_paired_reads(samfile):
    paired_reads = OrderedDict()
    with pysam.AlignmentFile(samfile, 'rb') as samfile:
        for idx, read in enumerate(samfile.fetch()):
            name = read.query_name
            paired_reads.setdefault(name, []).append(read)
    for header, pair in paired_reads.items():
        # if len(pair) > 2:
        #     raise RuntimeError(
        #         'Malformed SAM file: too many reads in one pair: {}'
        #         .format(header)
        #     )
        yield header, pair
