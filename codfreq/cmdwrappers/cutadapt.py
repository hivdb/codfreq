import os
import json
from subprocess import Popen, PIPE
from typing import List, Optional, TypedDict

from .base import raise_on_proc_error


class CutadaptConfig(TypedDict, total=False):
    adapter3: Optional[str]
    adapter5: Optional[str]
    adapter53: Optional[str]
    error_rate: Optional[float]
    no_indels: Optional[bool]
    times: Optional[int]
    min_overlap: Optional[int]


def load_config(
    config_path: str,
    adapter3_path: str,
    adapter5_path: str,
    adapter53_path: str
) -> Optional[CutadaptConfig]:
    if os.path.isfile(config_path):
        config: CutadaptConfig
        with open(config_path) as fp:
            config = json.load(fp)
        if os.path.isfile(adapter3_path):
            config['adapter3'] = 'file:' + adapter3_path
        if os.path.isfile(adapter5_path):
            config['adapter5'] = 'file:' + adapter5_path
        if os.path.isfile(adapter53_path):
            config['adapter53'] = 'file:' + adapter53_path
        return config
    else:
        return None


def cutadapt(
    merged_fastq_in: str,
    merged_fastq_out: str,
    adapter3: Optional[str] = None,
    adapter5: Optional[str] = None,
    adapter53: Optional[str] = None,
    error_rate: Optional[float] = None,
    no_indels: Optional[bool] = None,
    times: Optional[int] = None,
    min_overlap: Optional[int] = None
) -> None:
    command: List[str] = ['cutadapt', '-j', '0']
    if adapter3 is not None:
        command.extend(['-a', adapter3])
    if adapter5 is not None:
        command.extend(['-g', adapter5])
    if adapter53 is not None:
        command.extend(['-b', adapter53])
    if error_rate is not None:
        command.extend(['-e', str(error_rate)])
    if no_indels is not None:
        command.append('--no-indels')
    if times is not None:
        command.extend(['-n', str(times)])
    if min_overlap is not None:
        command.extend(['-O', str(min_overlap)])

    command.extend([
        '-o', merged_fastq_out,
        merged_fastq_in
    ])

    proc: Popen = Popen(
        command,
        stdout=PIPE,
        stderr=PIPE,
        encoding='U8')

    out: str
    error: str
    out, error = proc.communicate()

    raise_on_proc_error(proc, error)

    with open(merged_fastq_out + '.cutadapt.log', 'w') as fp:
        fp.write(' '.join(command) + '\n')
        fp.write(out)
        fp.write(error)
