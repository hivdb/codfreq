from typing import TextIO, List, Optional
from .codfreq_types import Sequence


def load(fp: TextIO) -> List[Sequence]:
    sequences: List[Sequence] = []
    header: Optional[str] = None
    curseq: str = ''
    for line in fp:
        if line.startswith('>'):
            if header and curseq:
                sequences.append({
                    'header': header,
                    'sequence': curseq.upper()
                })
            header = line[1:].strip()
            curseq = ''
        elif line.startswith('#'):
            continue
        else:
            curseq += line.strip()
    if header and curseq:
        sequences.append({
            'header': header,
            'sequence': curseq.upper()
        })
    return sequences
