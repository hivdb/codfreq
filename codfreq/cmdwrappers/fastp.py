import os
import json
from subprocess import Popen, PIPE
from typing import List, Optional, TypedDict
import multiprocessing

THREADS = '{}'.format(multiprocessing.cpu_count() // 2 + 1)


class FASTPConfig(TypedDict, total=False):
    disable_quality_filtering: Optional[bool]
    qualified_quality_phred: Optional[int]
    unqualified_percent_limit: Optional[int]
    average_qual: Optional[int]
    disable_length_filtering: Optional[bool]
    length_required: Optional[int]
    length_limit: Optional[int]
    disable_adapter_trimming: Optional[bool]
    adapter_sequence: Optional[str]
    adapter_sequence_r2: Optional[str]


def load_fastp_config(config_path: str) -> Optional[FASTPConfig]:
    if os.path.isfile(config_path):
        with open(config_path) as fp:
            config: FASTPConfig = json.load(fp)
            return config
    else:
        return None


def fastp(
    fastq1in: str,
    fastq2in: Optional[str],
    fastq1out: str,
    fastq2out: Optional[str],
    disable_quality_filtering: Optional[bool] = False,
    qualified_quality_phred: Optional[int] = None,
    unqualified_percent_limit: Optional[int] = None,
    average_qual: Optional[int] = None,
    disable_length_filtering: Optional[bool] = False,
    length_required: Optional[int] = None,
    length_limit: Optional[int] = None,
    disable_adapter_trimming: Optional[bool] = False,
    adapter_sequence: Optional[str] = None,
    adapter_sequence_r2: Optional[str] = None
) -> None:
    command: List[str] = [
        'fastp',
        '-w', str(THREADS),
        '-i', fastq1in,
        '-o', fastq1out
    ]
    if fastq2in and fastq2out:
        command.extend([
            '-I', fastq2in,
            '-O', fastq2out
        ])
    if disable_quality_filtering:
        command.append('-Q')
    else:
        if qualified_quality_phred is not None:
            command.extend(['-q', str(qualified_quality_phred)])
        if unqualified_percent_limit is not None:
            command.extend(['-u', str(unqualified_percent_limit)])
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

    if error:
        raise RuntimeError(error)

    with open(os.path.splitext(fastq1out)[0] + '.log', 'w') as fp:
        fp.write(out)
