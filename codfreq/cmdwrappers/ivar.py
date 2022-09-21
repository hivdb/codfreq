import os
import json
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile
from typing import List, Optional, TypedDict
import multiprocessing

from .base import execute, raise_on_proc_error

THREADS = '{}'.format(multiprocessing.cpu_count() // 2 + 1)


class TrimConfig(TypedDict, total=False):
    primers_bed: Optional[str]
    min_length: Optional[int]
    min_quality: Optional[int]
    sliding_window_width: Optional[int]
    include_reads_with_no_primers: Optional[bool]


def load_trim_config(
    config_path: str,
    primers_bed: str
) -> Optional[TrimConfig]:
    if os.path.isfile(config_path):
        config: TrimConfig
        with open(config_path) as fp:
            config = json.load(fp)
        if os.path.isfile(primers_bed):
            config['primers_bed'] = primers_bed
        return config
    else:
        return None


def trim(
    input_bam: str,
    output_bam: str,
    primers_bed: Optional[str] = None,
    min_length: Optional[int] = None,
    min_quality: Optional[int] = None,
    sliding_window_width: Optional[int] = None,
    include_reads_with_no_primers: Optional[bool] = False
) -> None:
    basedir, filename = os.path.split(input_bam)

    command: List[str] = [
        'ivar', 'trim',
        '-i', input_bam
    ]
    if primers_bed is not None:
        command.extend(['-b', primers_bed])
    if min_length is not None:
        command.extend(['-m', str(min_length)])
    if min_quality is not None:
        command.extend(['-q', str(min_quality)])
    if sliding_window_width is not None:
        command.extend(['-s', str(sliding_window_width)])
    if include_reads_with_no_primers:
        command.append('-e')

    with NamedTemporaryFile(suffix='.bam') as tmpbam:
        command.extend(['-p', os.path.splitext(tmpbam.name)[0]])

        proc: Popen = Popen(
            command,
            stdout=PIPE,
            stderr=PIPE,
            encoding='U8')

        out: str
        error: str
        out, error = proc.communicate()

        raise_on_proc_error(proc, error)

        out_samsort, err_samsort = execute([
            'samtools',
            'sort',
            '-@', THREADS,
            '-O', 'bam',
            '-o', output_bam,
            tmpbam.name
        ])
    out_samidx, err_samidx = execute([
        'samtools',
        'index',
        '-@', THREADS,
        output_bam
    ])

    with open(output_bam + '.ivar.log', 'w') as fp:
        fp.write(out)
        fp.write(error)
        fp.write(out_samsort)
        fp.write(err_samsort)
        fp.write(out_samidx)
        fp.write(err_samidx)
