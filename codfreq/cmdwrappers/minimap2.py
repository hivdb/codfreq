#! /usr/bin/env python

import os
from subprocess import Popen, PIPE

from .base import execute, raise_on_proc_error, refinit_func, align_func

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
def minimap2_align(refseq, fastq1, fastq2, bam):
    command = [
        'minimap2',
        *MINIMAP2_ARGS,
        '-a', refseq, fastq1
    ]
    if fastq2:
        command.append(fastq2)
    proc_minimap2 = Popen(
        command,
        stdout=PIPE,
        stderr=PIPE,
        encoding='U8')
    proc_sam2bam = Popen(
        ['samtools', 'sort', '-O', 'bam', '-o', bam],
        stdin=proc_minimap2.stdout,
        stdout=PIPE,
        stderr=PIPE,
        encoding='U8')
    proc_minimap2.stdout.close()
    err_minimap2 = proc_minimap2.stderr.read()
    proc_minimap2.stderr.close()
    raise_on_proc_error(proc_minimap2, err_minimap2)
    out_sam2bam, err_sam2bam = proc_sam2bam.communicate()
    raise_on_proc_error(proc_sam2bam, err_sam2bam)
    out_samidx, err_samidx = execute(['samtools', 'index', bam])
    with open(os.path.splitext(bam)[0] + '.log', 'w') as fp:
        fp.write(err_minimap2)
        fp.write(out_sam2bam)
        fp.write(err_sam2bam)
        fp.write(out_samidx)
        fp.write(err_samidx)
    return {'overall_rate': -1.}
