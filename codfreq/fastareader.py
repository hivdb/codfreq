from typing import TextIO
from .codfreq_types import Sequence


def load(fp: TextIO) -> list[Sequence]:
    """Load FASTA entries from a file-like object.

    :param fp: Input stream positioned at a FASTA file.
    :type fp: TextIO
    :returns: Parsed sequence objects.
    :rtype: list[Sequence]
    """

    sequences: list[Sequence] = []
    header: str | None = None
    curseq: bytearray = bytearray()
    for line in fp:
        if line.startswith('>'):
            if header and curseq:
                sequences.append(
                    Sequence(
                        header=header,
                        sequence=curseq.upper().decode('U8')
                    )
                )
            header = line[1:].strip()
            curseq = bytearray()
        elif line.startswith('#'):
            continue
        else:
            curseq.extend(
                line.strip().encode('ASCII', errors='ignore')
            )
    if header and curseq:
        sequences.append(
            Sequence(
                header=header,
                sequence=curseq.upper().decode('U8')
            )
        )
    return sequences
