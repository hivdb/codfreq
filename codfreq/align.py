#! /usr/bin/env python

import os
import re
import csv
import json
import click  # type: ignore
import tempfile
from itertools import combinations
from collections import defaultdict
from typing import (
    TextIO,
    Generator,
    List,
    Tuple,
    DefaultDict,
    Set
)

from .codfreq_types import Profile, PairedFASTQ, CodFreqRow

from .sam2codfreq import (
    sam2codfreq_all,
    CODFREQ_HEADER
)
from .sam2consensus import create_untrans_region_consensus
from .cmdwrappers import (
    get_programs, get_refinit, get_align
)
from .filename_helper import (
    suggest_pair_name,
    name_bamfile,
    name_codfreq
)
# from .reference_helper import get_refaas

ENCODING = 'UTF-8'
FILENAME_DELIMITERS = (' ', '_', '-')
PAIRED_FASTQ_MARKER = ('1', '2')
INVALID_PAIRED_FASTQ_MARKER = re.compile(r'01|02|10|20')


def find_paired_marker(text1: str, text2: str) -> int:
    pos: int
    a: str
    b: str
    diffcount: int = 0
    diffpos: int = -1
    if (
        INVALID_PAIRED_FASTQ_MARKER.search(text1) or
        INVALID_PAIRED_FASTQ_MARKER.search(text2)
    ):
        return -1

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


def find_paired_fastq_patterns(
    filenames: List[str]
) -> Generator[PairedFASTQ, None, None]:
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
    left: str
    right: str
    invalid: bool
    delimiter: str
    diffcount: int
    diffoffset: int
    reverse: int
    pattern: Tuple[
        str,  # delimiter
        int,  # diffoffset
        int,  # pos_paired_marker
        int,  # reverse
    ]
    pairs: List[Tuple[str, str]]
    patterns: DefaultDict[
        Tuple[
            str,  # delimiter
            int,  # diffoffset
            int,  # pos_paired_marker
            int,  # reverse
        ],
        List[Tuple[str, str]]
    ] = defaultdict(list)
    for fn1, fn2 in combinations(filenames, 2):
        if len(fn1) != len(fn2):
            continue
        for delimiter in FILENAME_DELIMITERS:
            if delimiter not in fn1 or delimiter not in fn2:
                continue
            chunks1: List[str] = fn1.split(delimiter)
            chunks2: List[str] = fn2.split(delimiter)
            if len(chunks1) != len(chunks2):
                continue
            for reverse in range(2):
                diffcount = 0
                diffoffset = -1
                invalid = False
                if reverse:
                    chunks1.reverse()
                    chunks2.reverse()
                for n, (left, right) in enumerate(zip(chunks1, chunks2)):
                    if diffcount > 1:
                        invalid = True
                        break
                    if left == right:
                        continue
                    pos_paired_marker: int = find_paired_marker(left, right)
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
    covered: Set[str] = set()
    for pattern, pairs in sorted(
            patterns.items(), key=lambda p: (-len(p[1]), -p[0][3])):
        known: Set[str] = set()
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
            for pair in pairs:
                yield {
                    'name': suggest_pair_name(pair, pattern),
                    'pair': pair,
                    'n': 2
                }
    if len(filenames) > len(covered):
        remains: List[str] = sorted(set(filenames) - covered)
        pattern = ('', -1, -1, -1)
        for left in remains:
            yield {
                'name': suggest_pair_name((left, None), pattern),
                'pair': (left, None),
                'n': 1
            }


def complete_paired_fastqs(
    paired_fastqs: Generator[PairedFASTQ, None, None],
    dirpath: str
) -> Generator[PairedFASTQ, None, None]:
    for pairobj in paired_fastqs:
        yield {
            'name': os.path.join(dirpath, pairobj['name']),
            'pair': tuple(
                os.path.join(dirpath, fn) if fn else None
                for fn in pairobj['pair']
            ),
            'n': pairobj['n']
        }


def find_paired_fastqs(workdir: str) -> Generator[PairedFASTQ, None, None]:
    pairinfo: str = os.path.join(workdir, 'pairinfo.json')
    if os.path.isfile(pairinfo):
        with open(pairinfo) as fp:
            yield from complete_paired_fastqs(
                json.load(fp),
                workdir
            )
    else:
        for dirpath, _, filenames in os.walk(workdir, followlinks=True):
            filenames = [
                fn for fn in filenames
                if fn[-6:].lower() == '.fastq'
                or fn[-9:].lower() == '.fastq.gz'
            ]
            yield from complete_paired_fastqs(
                find_paired_fastq_patterns(filenames),
                dirpath
            )


def align_with_profile(
    paired_fastqs: List[PairedFASTQ],
    program: str,
    profile: Profile,
    log_format: str
) -> None:
    with tempfile.TemporaryDirectory('codfreq') as tmpdir:
        refpath = os.path.join(tmpdir, 'ref.fas')
        refinit = get_refinit(program)
        alignfunc = get_align(program)
        for config in profile['fragmentConfig']:
            if 'refSequence' not in config:
                continue
            refname = config['fragmentName']
            refseq = config['refSequence']
            with open(refpath, 'w') as fp:
                fp.write('>{}\n{}\n\n'.format(refname, refseq))
        # refaas = get_refaas(reference)
            for pairobj in paired_fastqs:
                samfile = name_bamfile(pairobj['name'], refname)
                refinit(refpath)
                if log_format == 'text':
                    click.echo(
                        'Aligning {} with {}...'
                        .format(pairobj['name'], refname)
                    )
                else:
                    click.echo(json.dumps({
                        'op': 'alignment',
                        'status': 'working',
                        'query': pairobj['name'],
                        'target': refname
                    }))
                alignfunc(refpath, *pairobj['pair'], samfile)
                if log_format == 'text':
                    click.echo('Done')
                else:
                    click.echo(json.dumps({
                        'op': 'alignment',
                        'status': 'done',
                        'query': pairobj['name'],
                        'target': refname
                    }))


def align(
    workdir: str,
    program: str,
    profile: TextIO,
    log_format: str
) -> None:
    row: CodFreqRow
    profile_obj: Profile = json.load(profile)
    paired_fastqs = list(find_paired_fastqs(workdir))
    align_with_profile(paired_fastqs, program, profile_obj, log_format)

    for pairobj in paired_fastqs:
        codfreqfile = name_codfreq(pairobj['name'])
        with open(codfreqfile, 'w', encoding='utf-8-sig') as fp:
            writer = csv.DictWriter(fp, CODFREQ_HEADER)
            writer.writeheader()
            for row in sam2codfreq_all(
                name=pairobj['name'],
                fnpair=pairobj['pair'],
                profile=profile_obj,
                log_format=log_format
            ):
                writer.writerow({
                    **row,
                    'codon': row['codon'].decode(ENCODING)
                })

        create_untrans_region_consensus(pairobj['name'], profile_obj)


@click.command()
@click.argument(
    'workdir',
    type=click.Path(exists=True, file_okay=False,
                    dir_okay=True, resolve_path=True))
@click.option(
    '-p', '--program',
    required=True,
    type=click.Choice(get_programs()))
@click.option(
    '-r', '--profile',
    required=True,
    type=click.File('r', encoding=ENCODING))
@click.option(
    '--log-format',
    type=click.Choice(['text', 'json']),
    default='text', show_default=True)
@click.option(
    '--enable-profiling/--disable-profiling',
    default=False,
    help='Enable cProfile')
def align_cmd(
    workdir: str,
    program: str,
    profile: TextIO,
    log_format: str,
    enable_profiling: bool
) -> None:
    if enable_profiling:
        import cProfile
        import pstats
        profile_obj = None
        try:
            with cProfile.Profile() as profile_obj:
                align(workdir, program, profile, log_format)
        finally:
            if profile_obj is not None:
                ps = pstats.Stats(profile_obj)
                ps.print_stats()
    else:
        align(workdir, program, profile, log_format)


if __name__ == '__main__':
    align_cmd()
