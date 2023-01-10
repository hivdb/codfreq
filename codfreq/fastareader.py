from collections import defaultdict
from typing import TextIO, List, Optional, Tuple, Dict

from .posnas import PosNA
from .codfreq_types import Sequence

MISSING: int = ord('-')


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


def write_multi_alignment(
    fh: TextIO,
    seqs: List[Tuple[str, Tuple[Optional[PosNA], ...]]]
) -> None:
    pos_sizes: Dict[int, int] = defaultdict(int)
    for header, seq in seqs:
        for node in seq:
            if node is None:
                continue
            if node.bp >= pos_sizes[node.pos]:
                pos_sizes[node.pos] = node.bp + 1
    for header, seq in seqs:
        fh.write(f'>{header}\n')
        seqmap: Dict[int, List[int]] = {
            pos: [MISSING] * size for pos, size in pos_sizes.items()
        }
        for node in seq:
            if node is None:
                continue
            seqmap[node.pos][node.bp] = node.na

        for _, nas in sorted(seqmap.items()):
            fh.write(bytes(nas).decode('ASCII'))
        fh.write('\n')
