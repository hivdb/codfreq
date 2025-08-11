"""Alignment utilities for preparing CodFreq inputs."""

import os
import re
import csv
import json
import sys
import tempfile
import multiprocessing
from itertools import combinations
from collections import defaultdict
from pathlib import Path
from typing import (
    TextIO,
    Annotated,
)
from collections.abc import Iterable, Generator

import rich
import typer

from .codfreq_types import Profile, PairedFASTQ, CodFreqRow

from .sam2codfreq import (
    sam2codfreq_all,
    CODFREQ_HEADER
)
from .sam2consensus import create_untrans_region_consensus
from .cmdwrappers import (
    fastp, cutadapt, ivar, get_refinit, get_align
)
from .enums import LogFormat, Program
from .filename_helper import (
    suggest_pair_name,
    name_bamfile,
    name_codfreq
)

ENCODING = 'UTF-8'
REQUIRED_PROFILE_VERSION = '20221213'
FILENAME_DELIMITERS = (' ', '_', '-')
PAIRED_FASTQ_MARKER = ('1', '2')
INVALID_PAIRED_FASTQ_MARKER = re.compile(r'[1-9]0*[12]|[^0]00+[12]|[12]\d')

app = typer.Typer(pretty_exceptions_enable=False)


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
    filenames: list[str],
    autopairing: bool
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
    pattern: tuple[
        str,  # delimiter
        int,  # diffoffset
        int,  # pos_paired_marker
        int,  # reverse
    ]
    pairs: list[tuple[str, str]]
    patterns: defaultdict[
        tuple[
            str,  # delimiter
            int,  # diffoffset
            int,  # pos_paired_marker
            int,  # reverse
        ],
        list[tuple[str, str]]
    ] = defaultdict(list)
    if autopairing:
        for fn1, fn2 in combinations(filenames, 2):
            if len(fn1) != len(fn2):
                continue
            for delimiter in FILENAME_DELIMITERS:
                if delimiter not in fn1 or delimiter not in fn2:
                    continue
                chunks1: list[str] = fn1.split(delimiter)
                chunks2: list[str] = fn2.split(delimiter)
                if len(chunks1) != len(chunks2):
                    continue  # pragma: no cover - unequal chunk counts
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
                        pos_paired_marker: int = \
                            find_paired_marker(left, right)
                        if pos_paired_marker < 0:
                            invalid = True
                            break
                        diffoffset = n
                        diffcount += 1
                    if not invalid:
                        if fn1 > fn2:
                            # sort by filename
                            fn1, fn2 = fn2, fn1  # pragma: no cover - swap
                        patterns[(
                            delimiter,
                            diffoffset,
                            pos_paired_marker,
                            reverse
                        )].append((fn1, fn2))
    covered: set[str] = set()
    if autopairing:
        for pattern, pairs in sorted(
                patterns.items(), key=lambda p: (-len(p[1]), -p[0][3])):
            known: set[str] = set()
            invalid = False
            for left, right in pairs:
                if left in covered or right in covered:
                    # a pattern is invalid if the pairs is already matched
                    # by a previous pattern
                    invalid = True
                    break

                if left in known or right in known:
                    # a pattern is invalid if there's duplicate in pairs
                    invalid = True  # pragma: no cover - repeated chunks
                    break  # pragma: no cover - stop after duplicate
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
        remains: list[str] = sorted(set(filenames) - covered)
        pattern = ('', -1, -1, -1)
        for left in remains:
            yield {
                'name': suggest_pair_name((left, None), pattern),
                'pair': (left, None),
                'n': 1
            }


def complete_paired_fastqs(
    paired_fastqs: Iterable[PairedFASTQ],
    dirpath: str
) -> Generator[PairedFASTQ, None, None]:
    for pairobj in paired_fastqs:
        yield {
            'name': os.path.join(dirpath, pairobj['name']),
            'pair': (
                os.path.join(dirpath, pairobj['pair'][0]),
                os.path.join(dirpath, pairobj['pair'][1])
                if pairobj['pair'][1] else None
            ),
            'n': pairobj['n']
        }


def find_paired_fastqs(
    workdir: str,
    autopairing: bool
) -> Generator[PairedFASTQ, None, None]:
    pairinfo: str = os.path.join(workdir, 'pairinfo.json')
    if os.path.isfile(pairinfo):
        with open(pairinfo) as fp:
            yield from complete_paired_fastqs(
                json.load(fp),
                workdir
            )
    else:
        pairinfo_list: list[PairedFASTQ] = []
        for dirpath, _, filenames in os.walk(workdir, followlinks=True):
            filenames = [
                fn for fn in filenames
                if (
                    fn[-6:].lower() == '.fastq'
                    or fn[-9:].lower() == '.fastq.gz'
                ) and not (
                    fn[-12:].lower() == 'merged.fastq'
                    or fn[-15:].lower() == 'merged.fastq.gz'
                )
            ]
            rel_dirpath = os.path.relpath(dirpath, workdir)
            pairinfo_list.extend(complete_paired_fastqs(
                find_paired_fastq_patterns(filenames, autopairing),
                rel_dirpath
            ))
        with open(pairinfo, 'w') as fp:
            json.dump(pairinfo_list, fp, indent=2)
        yield from complete_paired_fastqs(pairinfo_list, workdir)


def fastp_preprocess(
    paired_fastq: PairedFASTQ,
    fastp_config: fastp.FASTPConfig,
    log_format: LogFormat
) -> PairedFASTQ:
    """Merge paired FASTQ reads using fastp.

    :param paired_fastq: Input paired FASTQ metadata.
    :param fastp_config: Configuration options for fastp.
    :param log_format: Output format for logging.
    :type log_format: LogFormat
    :returns: Metadata for the merged FASTQ file.
    """
    if log_format == LogFormat.text:
        rich.print(
            'Pre-processing {} using fastp...'
            .format(paired_fastq['name'])
        )
    else:
        print(json.dumps({
            'op': 'preprocess',
            'status': 'working',
            'query': paired_fastq['name']
        }))
    merged_fastq: PairedFASTQ = {
        'name': paired_fastq['name'],
        'pair': (
            os.path.join(
                os.path.dirname(paired_fastq['pair'][0]),
                '{}.merged.fastq.gz'.format(paired_fastq['name'])
            ),
            None
        ),
        'n': 1
    }
    fastp.fastp(
        paired_fastq['pair'][0],
        paired_fastq['pair'][1],
        merged_fastq['pair'][0],
        **fastp_config
    )
    if log_format == LogFormat.text:
        rich.print('Done')
    else:
        print(json.dumps({
            'op': 'preprocess',
            'status': 'done',
            'query': paired_fastq['name']
        }))
    return merged_fastq


def ivar_trim(
    input_bam: str,
    output_bam: str,
    ivar_trim_config: ivar.TrimConfig,
    log_format: LogFormat
) -> None:
    """Trim primer sequences using ivar.

    :param input_bam: Input BAM file path.
    :param output_bam: Output trimmed BAM file path.
    :param ivar_trim_config: Configuration for ivar trim.
    :param log_format: Logging output format.
    :type log_format: LogFormat
    :returns: None
    """
    name: str = os.path.basename(input_bam)
    if log_format == LogFormat.text:
        rich.print(
            'Trimming {} using ivar...'
            .format(name)
        )
    else:
        print(json.dumps({
            'op': 'trim',
            'status': 'working',
            'command': 'ivar',
            'query': name
        }))
    ivar.trim(
        input_bam,
        output_bam,
        **ivar_trim_config
    )
    if log_format == LogFormat.text:
        rich.print('Done')
    else:
        print(json.dumps({
            'op': 'trim',
            'status': 'done',
            'command': 'ivar',
            'query': name
        }))


def cutadapt_trim(
    merged_fastq: PairedFASTQ,
    cutadapt_config: cutadapt.CutadaptConfig,
    log_format: LogFormat
) -> PairedFASTQ:
    """Trim reads using cutadapt.

    :param merged_fastq: FASTQ metadata after merging.
    :param cutadapt_config: Configuration for cutadapt.
    :param log_format: Logging output format.
    :type log_format: LogFormat
    :returns: Metadata for the trimmed FASTQ file.
    """
    name: str = merged_fastq['name']
    output_fastq: PairedFASTQ = {
        'name': merged_fastq['name'],
        'pair': (
            os.path.join(
                os.path.dirname(merged_fastq['pair'][0]),
                '{}.merged-trimed.fastq.gz'.format(merged_fastq['name'])
            ),
            None
        ),
        'n': 1
    }
    if log_format == LogFormat.text:
        rich.print(
            'Trimming {} using cutadapt...'
            .format(name)
        )
    else:
        print(json.dumps({
            'op': 'trim',
            'status': 'working',
            'command': 'cutadapt',
            'query': name
        }))
    cutadapt.cutadapt(
        merged_fastq['pair'][0],
        output_fastq['pair'][0],
        **cutadapt_config
    )
    if log_format == LogFormat.text:
        rich.print('Done')
    else:
        print(json.dumps({
            'op': 'trim',
            'status': 'done',
            'command': 'cutadapt',
            'query': name
        }))
    return output_fastq


def align_with_profile(
    paired_fastq: PairedFASTQ,
    program: Program,
    profile: Profile,
    log_format: LogFormat,
    fastp_config: fastp.FASTPConfig,
    cutadapt_config: cutadapt.CutadaptConfig | None,
    ivar_trim_config: ivar.TrimConfig | None
) -> None:
    """Align reads to references defined in the profile.

    :param paired_fastq: FASTQ pair to align.
    :param program: Alignment program to execute.
    :type program: Program
    :param profile: Profile configuration object.
    :param log_format: Logging output format.
    :type log_format: LogFormat
    :param fastp_config: ``fastp`` preprocessing settings.
    :param cutadapt_config: ``cutadapt`` trimming settings, if any.
    :param ivar_trim_config: ``ivar`` trimming settings, if any.
    :returns: None
    """
    paired_fastq = fastp_preprocess(paired_fastq, fastp_config, log_format)

    if cutadapt_config is not None:
        paired_fastq = cutadapt_trim(
            paired_fastq, cutadapt_config, log_format)

    with tempfile.TemporaryDirectory('codfreq') as tmpdir:
        refpath = os.path.join(tmpdir, 'ref.fas')
        refinit = get_refinit(program.value)
        alignfunc = get_align(program.value)
        for config in profile.fragmentConfig:
            if config.refSequence is None:
                continue  # pragma: no cover - missing refSequence
            refname = config.fragmentName
            refseq = config.refSequence
            with open(refpath, 'w') as fp:
                fp.write(f'>{refname}\n{refseq}\n\n')

            orig_bamfile = name_bamfile(
                paired_fastq['name'],
                refname,
                is_trimmed=False)
            trimmed_bamfile = name_bamfile(
                paired_fastq['name'],
                refname,
                is_trimmed=True)
            refinit(refpath)
            if log_format == LogFormat.text:
                rich.print(
                    'Aligning {} with {}...'
                    .format(paired_fastq['name'], refname)
                )
            else:
                print(json.dumps({
                    'op': 'alignment',
                    'status': 'working',
                    'query': paired_fastq['name'],
                    'target': refname
                }))
            alignfunc(refpath, *paired_fastq['pair'], orig_bamfile)
            if log_format == LogFormat.text:
                rich.print('Done')
            else:
                print(json.dumps({
                    'op': 'alignment',
                    'status': 'done',
                    'query': paired_fastq['name'],
                    'target': refname
                }))
            if ivar_trim_config is None:
                os.replace(orig_bamfile, trimmed_bamfile)
                os.replace(orig_bamfile + '.bai', trimmed_bamfile + '.bai')
                os.replace(orig_bamfile + '.minimap2.log',
                           trimmed_bamfile + '.minimap2.log')
            else:
                ivar_trim(
                    orig_bamfile,
                    trimmed_bamfile,
                    ivar_trim_config,
                    log_format)


def align(
    workdir: Path,
    program: Program,
    profile: TextIO,
    workers: int,
    log_format: LogFormat,
    autopairing: bool
) -> None:
    """Run the alignment pipeline to produce CodFreq files.

    :param workdir: Working directory containing inputs and outputs.
    :type workdir: Path
    :param program: Alignment program to execute.
    :type program: Program
    :param profile: Open profile configuration file.
    :param workers: Number of worker processes to use.
    :param log_format: Logging output format.
    :type log_format: LogFormat
    :param autopairing: Enable automatic FASTQ pairing.
    :returns: None
    :raises typer.Abort: If the profile version is incompatible.
    """
    row: CodFreqRow
    profile_data = json.load(profile)
    if profile_data.get('version') != REQUIRED_PROFILE_VERSION:
        rich.print(
            'Incompatible profile detected. Download the latest profile files '
            'from: https://github.com/hivdb/codfreq/tree/main/profiles',
            file=sys.stderr)
        raise typer.Abort()
    profile_obj = Profile.model_validate(profile_data)
    paired_fastqs = list(find_paired_fastqs(str(workdir), autopairing))

    fastp_config: fastp.FASTPConfig = fastp.load_config(
        str(workdir / 'fastp-config.json')
    )
    cutadapt_config: cutadapt.CutadaptConfig | None = cutadapt.load_config(
        str(workdir / 'cutadapt-config.json'),
        adapter3_path=str(workdir / 'primers3.fa'),
        adapter5_path=str(workdir / 'primers5.fa'),
        adapter53_path=str(workdir / 'primers53.fa')
    )
    ivar_trim_config: ivar.TrimConfig | None = ivar.load_trim_config(
        str(workdir / 'ivar-trim-config.json'),
        str(workdir / 'primers.bed')
    )

    for pairobj in paired_fastqs:
        align_with_profile(
            pairobj,
            program,
            profile_obj,
            log_format,
            fastp_config=fastp_config,
            cutadapt_config=cutadapt_config,
            ivar_trim_config=ivar_trim_config
        )
        codfreqfile = name_codfreq(pairobj['name'])
        with open(codfreqfile, 'w', encoding='utf-8-sig') as fp:
            writer = csv.DictWriter(fp, CODFREQ_HEADER)
            writer.writeheader()
            for row in sam2codfreq_all(
                name=pairobj['name'],
                fnpair=pairobj['pair'],
                profile=profile_obj,
                workers=workers,
                log_format=log_format
            ):
                writer.writerow({
                    **row,
                    'codon': row['codon'].decode(ENCODING)
                })

        create_untrans_region_consensus(
            pairobj['name'],
            profile_obj
        )


@app.command()
def align_cmd(
    workdir: Annotated[
        Path,
        typer.Argument(
            ..., exists=True, file_okay=False, dir_okay=True, resolve_path=True
        ),
    ],
    program: Annotated[
        Program,
        typer.Option('-p', help='Alignment program'),
    ],
    profile: Annotated[
        typer.FileText,
        typer.Option('-r', encoding=ENCODING, help='Profile JSON file'),
    ],
    log_format: Annotated[
        LogFormat,
        typer.Option(
            '--log-format',
            show_default=True,
            help='Log output format',
        ),
    ] = LogFormat.text,
    enable_profiling: Annotated[
        bool,
        typer.Option(
            '--enable-profiling/--disable-profiling',
            help='Enable/disable cProfile',
        ),
    ] = False,
    autopairing: Annotated[
        bool,
        typer.Option(
            '--autopairing/--no-autopairing',
            help='Enable/disable automatical FASTQ pairing algorithm',
        ),
    ] = True,
    workers: Annotated[
        int,
        typer.Option(
            '--workers',
            show_default=True,
            help='Number of sub-process workers to be used',
        ),
    ] = multiprocessing.cpu_count(),
) -> None:
    """Command-line interface for :func:`align`.

    :param workdir: Working directory with FASTQ files.
    :param program: Alignment program to execute.
    :type program: Program
    :param profile: Profile configuration file handle.
    :param log_format: Format for log output.
    :type log_format: LogFormat
    :param enable_profiling: Run with cProfile if ``True``.
    :param autopairing: Automatically pair FASTQ files if ``True``.
    :param workers: Number of worker processes.
    :returns: None
    """
    if enable_profiling:
        import cProfile
        import pstats
        profile_obj = None
        try:
            with cProfile.Profile() as profile_obj:
                align(
                    workdir,
                    program,
                    profile,
                    workers,
                    log_format,
                    autopairing,
                )
        finally:
            if profile_obj is not None:
                ps = pstats.Stats(profile_obj)
                ps.print_stats()
    else:
        align(
            workdir,
            program,
            profile,
            workers,
            log_format,
            autopairing,
        )


if __name__ == '__main__':
    app()  # pragma: no cover - manual CLI execution
