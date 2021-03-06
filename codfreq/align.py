#! /usr/bin/env python

# from . import rayinit  # noqa
import os
import csv
import click
from itertools import combinations
from collections import defaultdict

from .sam2codfreq import sam2codfreq, CODFREQ_HEADER
from .cmdwrappers import (
    get_programs, get_refinit, get_align
)
from .filename_helper import (
    name_samfile,
    replace_ext
)
from .reference_helper import get_refaas

FILENAME_DELIMITERS = (' ', '_', '-')
PAIRED_FASTQ_MARKER = ('1', '2')


def find_paired_marker(text1, text2):
    diffcount = 0
    diffpos = -1
    for pos, (a, b) in enumerate(zip(text1, text2)):
        if diffcount > 1:
            return -1
        if a == b:
            continue
        if a not in PAIRED_FASTQ_MARKER or b not in PAIRED_FASTQ_MARKER:
            return -1
        diffcount += 1
        diffpos = pos
    return diffpos


def find_paired_fastq_patterns(filenames):
    """Smartly find paired FASTQ file patterns

    A valid filename pattern must meet:
    - use one of the valid delimiters (" ", "_" or "-") to separate the
      filename into different chunks
    - in one and only one chunk, a fixed position character changed from "1" to
      "2"

    Valid pair pattern examples:
      14258F_L001_R1_001.fastq.gz <-> 14258F_L001_R2_001.fastq.gz
      SampleExample_1.fastq <-> SampleExample_2.fastq

    Invalid pair pattern examples:
      SampleExample1.fastq <-> SampleExample2.fastq
      SampleExample_1.fastq <-> SampleExample_2.fastq.gz
      SampleExample_1.FASTQ.GZ <-> SampleExample_2.fastq.gz

    """
    patterns = defaultdict(list)
    for fn1, fn2 in combinations(filenames, 2):
        if len(fn1) != len(fn2):
            continue
        for delimiter in FILENAME_DELIMITERS:
            if delimiter not in fn1 or delimiter not in fn2:
                continue
            chunks1 = fn1.split(delimiter)
            chunks2 = fn2.split(delimiter)
            if len(chunks1) != len(chunks2):
                continue
            for reverse in range(2):
                diffcount = 0
                invalid = False
                diffoffset = -1
                if reverse:
                    chunks1.reverse()
                    chunks2.reverse()
                for n, (left, right) in enumerate(zip(chunks1, chunks2)):
                    if diffcount > 1:
                        invalid = True
                        break
                    if left == right:
                        continue
                    pos_paired_marker = find_paired_marker(left, right)
                    if pos_paired_marker < 0:
                        invalid = True
                        break
                    diffoffset = n
                    diffcount += 1
                if not invalid:
                    patterns[(
                        delimiter,
                        diffoffset,
                        pos_paired_marker,
                        reverse
                    )].append((fn1, fn2))
    covered = set()
    for pattern, pairs in sorted(
            patterns.items(), key=lambda p: (-len(p[1]), -p[0][3])):
        known = set()
        invalid = False
        for left, right in pairs:
            if left in covered or right in covered:
                # a pattern is invalid if the pairs is already matched
                # by a previous pattern
                invalid = True
                break

            if left in known or right in known:
                # a pattern is invalid if there's duplicate in pairs
                invalid = True
                break
            known.add(left)
            known.add(right)

        if not invalid:
            covered |= known
            yield pattern, pairs
    if len(filenames) > len(covered):
        remains = sorted(set(filenames) - covered)
        yield (None, -1, -1, -1), [(left, None) for left in remains]


def find_paired_fastqs(workdir):
    for dirpath, _, filenames in os.walk(workdir, followlinks=True):
        filenames = [
            fn for fn in filenames
            if fn[-6:].lower() == '.fastq'
            or fn[-9:].lower() == '.fastq.gz'
        ]
        yield from (
            (tuple(os.path.join(dirpath, fn) if fn else None
                   for fn in fnpair), pattern)
            for pattern, fnpairs in find_paired_fastq_patterns(filenames)
            for fnpair in fnpairs
        )


@click.command()
@click.argument(
    'workdir',
    type=click.Path(exists=True, file_okay=False,
                    dir_okay=True, resolve_path=True))
@click.option(
    '-g', '--gene',
    required=True,
    type=str,
    help='Gene name')
@click.option(
    '-p', '--program',
    required=True,
    type=click.Choice(get_programs()))
@click.option(
    '-r', '--reference',
    required=True,
    type=click.Path(exists=True, file_okay=True,
                    dir_okay=False, resolve_path=True))
@click.option(
    '--log-format',
    type=click.Choice(['text', 'json']),
    default='text', show_default=True)
def align(workdir, gene, program, reference, log_format):
    refinit = get_refinit(program)
    alignfunc = get_align(program)
    refaas = get_refaas(reference)
    for fnpair, pattern in find_paired_fastqs(workdir):
        samfile = name_samfile(fnpair, pattern, gene)
        refinit(reference)
        alignfunc(reference, *fnpair, samfile)
        codfreqfile = replace_ext(samfile, '.codfreq')
        with open(codfreqfile, 'w', encoding='utf-8-sig') as fp:
            writer = csv.DictWriter(fp, CODFREQ_HEADER)
            writer.writeheader()
            writer.writerows(sam2codfreq(
                samfile,
                gene=gene,
                refaas=refaas,
                reference_start=1,
                fnpair=fnpair,
                log_format=log_format
            ))


if __name__ == '__main__':
    align()
