#! /usr/bin/env python

import os
import re
from typing import List, Dict, Optional, Tuple

from .base import execute, refinit_func, align_func

BOWTIE2_ARGS: List[str] = [
    '--local',
    '--threads', '3',
    '--rdg', '15,3',  # read gap open, extend penalties
    '--rfg', '15,3',  # reference gap open, extend penalties
    '--no-unal'
    # '--ma', '2',     # match bonus
    # '--mp', '6',     # max mismatch penalty; lower qual = lower penalty
]


@refinit_func('bowtie2')
def bowtie2_refinit(refseq: str) -> None:
    refpath: str
    ext: str
    suffix: str
    refdir: str
    refname: str
    refpath, ext = os.path.splitext(refseq)
    suffixes: Tuple[str, ...] = (
        '1.bt2',
        '2.bt2',
        '3.bt2',
        '4.bt2',
        'rev.1.bt2',
        'rev.2.bt2'
    )
    if not all(os.path.isfile(
        os.path.extsep.join((refpath, suffix))
    ) for suffix in suffixes):
        # the reference index files are not previously generated
        refdir, refname = os.path.split(refpath)
        execute([
            'bowtie2-build',
            refpath + ext,
            refpath
        ])


@align_func('bowtie2')
def bowtie2_align(
    refseq: str,
    fastq1: str,
    fastq2: str,
    sam: str
) -> Dict[str, float]:
    refpath: str
    refdir: str
    refname: str
    logs: str
    refpath, _ = os.path.splitext(refseq)
    refdir, refname = os.path.split(refpath)
    command: List[str] = [
        'bowtie2',
        *BOWTIE2_ARGS,
        '-x', refpath,
        '-S', sam,
        '-U', fastq1
    ]
    if fastq2:
        command.extend(['-U', fastq2])
    logs, _ = execute(command)
    overall_rate: float = -1.
    overall_rate_match: Optional[re.Match] = re.search(
        r'\n(\d+\.\d+)% overall alignment rate',
        logs)
    if overall_rate_match:
        overall_rate_text = overall_rate_match.group(1)
        if overall_rate_text:
            overall_rate = float(overall_rate_text)
    with open(os.path.splitext(sam)[0] + '.log', 'w') as fp:
        fp.write(logs)
    return {'overall_rate': overall_rate}
