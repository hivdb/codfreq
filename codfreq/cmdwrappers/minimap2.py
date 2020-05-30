#! /usr/bin/env python

import os
from .base import execute, refinit_func, align_func

MINIMAP2_ARGS = [
    '-A', '2',        # matching score [2]
    '-B', '4',        # mismatch penalty [4]
    '-O', '4,24',     # gap open penalty [4,24]
    '-E', '2,1',      # gap extension penalty [2,1]
    '-z', '400,200',  # Z-drop score and inversion Z-drop score [400,200]
    '-s', '80',       # minimal peak DP alignment score [80]
    '-u', 'n',        # how to find GT-AG [n]
    '--sam-hit-only'
]


@refinit_func('minimap2')
def minimap2_refinit(refseq):
    return


@align_func('minimap2')
def minimap2_align(refseq, fastq1, fastq2, sam):
    command = [
        'minimap2',
        *MINIMAP2_ARGS,
        '-a', '-o', sam, refseq, fastq1
    ]
    if fastq2:
        command.append(fastq2)
    _, logs = execute(command)
    with open(os.path.splitext(sam)[0] + '.log', 'w') as fp:
        fp.write(logs)
    return {'overall_rate': -1.}
