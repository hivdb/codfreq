from typing import TextIO, List, Optional
from .codfreq_types import Sequence


def load(fp: TextIO) -> List[Sequence]:
    sequences: List[Sequence] = []
    header: Optional[str] = None
    curseq: bytearray = bytearray()
    for line in fp:
        if line.startswith('>'):
            if header and curseq:
                sequences.append({
                    'header': header,
                    'sequence': curseq.upper().decode('U8')
                })
            header = line[1:].strip()
            curseq = bytearray()
        elif line.startswith('#'):
            continue
        else:
            curseq.extend(
                line.strip()
                .encode('ASCII', errors='ignore')
            )
    if header and curseq:
        sequences.append({
            'header': header,
            'sequence': curseq.upper().decode('U8')
        })
    return sequences
