import os
import json
from subprocess import Popen, PIPE
from typing import List, Optional, TypedDict
import multiprocessing

from .base import raise_on_proc_error

THREADS = '{}'.format(multiprocessing.cpu_count() // 2 + 1)


class FASTPConfig(TypedDict, total=False):
    include_unmerged: Optional[bool]
    qualified_quality_phred: Optional[int]
    unqualified_percent_limit: Optional[int]
    n_base_limit: Optional[int]
    average_qual: Optional[int]
    length_required: Optional[int]
    length_limit: Optional[int]
    adapter_sequence: Optional[str]
    adapter_sequence_r2: Optional[str]
    disable_adapter_trimming: Optional[bool]
    disable_trim_poly_g: Optional[bool]
    disable_quality_filtering: Optional[bool]
    disable_length_filtering: Optional[bool]


def load_config(config_path: str) -> FASTPConfig:
    if os.path.isfile(config_path):
        with open(config_path) as fp:
            config: FASTPConfig = json.load(fp)
            return config
    else:
        return {}


def fastp(
    fastq1in: str,
    fastq2in: Optional[str],
    fastq_merged_out: str,
    include_unmerged: Optional[bool] = False,
    qualified_quality_phred: Optional[int] = None,
    unqualified_percent_limit: Optional[int] = None,
    n_base_limit: Optional[int] = None,
    average_qual: Optional[int] = None,
    length_required: Optional[int] = None,
    length_limit: Optional[int] = None,
    disable_adapter_trimming: Optional[bool] = False,
    disable_trim_poly_g: Optional[bool] = False,
    disable_quality_filtering: Optional[bool] = False,
    disable_length_filtering: Optional[bool] = False,
    adapter_sequence: Optional[str] = None,
    adapter_sequence_r2: Optional[str] = None
) -> None:
    command: List[str] = [
        'fastp',
        '-w', str(THREADS),
        '-i', fastq1in
    ]
    if fastq2in:
        command.extend([
            '-I', fastq2in,
            '-m', '--merged_out', fastq_merged_out,
        ])
        if include_unmerged:
            command.append('--include_unmerged')
    else:
        command.extend([
            '-o', fastq_merged_out
        ])
    if disable_quality_filtering:
        command.append('-Q')
    else:
        if qualified_quality_phred is not None:
            command.extend(['-q', str(qualified_quality_phred)])
        if unqualified_percent_limit is not None:
            command.extend(['-u', str(unqualified_percent_limit)])
        if n_base_limit is not None:
            command.extend(['-n', str(n_base_limit)])
        if average_qual is not None:
            command.extend(['-e', str(average_qual)])
    if disable_length_filtering:
        command.append('-L')
    else:
        if length_required is not None:
            command.extend(['-l', str(length_required)])
        if length_limit is not None:
            command.extend(['--length_limit', str(length_limit)])
    if disable_adapter_trimming:
        command.append('-A')
    else:
        if adapter_sequence:
            command.extend(['-a', adapter_sequence])
        if adapter_sequence_r2:
            command.extend(['--adapter_sequence_r2', adapter_sequence_r2])

    proc: Popen = Popen(
        command,
        stdout=PIPE,
        stderr=PIPE,
        encoding='U8')

    out: str
    error: str
    out, error = proc.communicate()

    raise_on_proc_error(proc, error)

    with open(os.path.splitext(fastq_merged_out)[0] + '.fastp.log', 'w') as fp:
        fp.write(' '.join(command) + '\n')
        fp.write(out)
        fp.write(error)
