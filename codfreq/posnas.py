import pysam  # type: ignore
import cython  # type: ignore
from tqdm import tqdm  # type: ignore
from array import array
from pysam import AlignedSegment  # type: ignore
from typing import List, Tuple, Optional, Generator

from .codfreq_types import NAPos, NAChar, SeqText, Header

ENCODING: str = 'UTF-8'
GAP: int = ord(b'-')

PosNA = Tuple[
    NAPos,   # refpos
    int,     # insertion_index
    NAChar,  # na
    int      # qua
]


@cython.ccall
@cython.inline
@cython.returns(list)
def iter_single_read_posnas(
    seq: SeqText,
    qua: array,
    aligned_pairs: List[Tuple[Optional[NAPos], Optional[NAPos]]]
) -> List[PosNA]:
    seqpos0: Optional[NAPos]
    refpos0: Optional[NAPos]
    refpos: NAPos
    insidx: int = 0
    n: NAChar
    q: int

    seqchars: bytes = bytes(seq, ENCODING)

    prev_refpos: int = 0
    prev_seqpos0: int = 0

    buffer_size: int = 0

    posnas: List[PosNA] = []

    for seqpos0, refpos0 in aligned_pairs:

        if refpos0 is None:
            # insertion
            refpos = prev_refpos
            insidx += 1
        else:
            refpos = refpos0 + 1
            insidx = 0
            prev_refpos = refpos

        if seqpos0 is None:
            # deletion
            n, q = GAP, qua[prev_seqpos0]
        else:
            n, q = seqchars[seqpos0], qua[seqpos0]
            prev_seqpos0 = seqpos0

        if refpos == 0:
            # insertion before the first ref position
            continue

        posnas.append((refpos, insidx, n, q))

        if insidx > 0:
            buffer_size += 1
        else:
            buffer_size = 0

    return posnas[:len(posnas) - buffer_size]


def iter_posnas(
    samfile: str,
    site_quality_cutoff: int = 0,
    progress: bool = True
) -> Generator[
    Tuple[str, List[PosNA]],
    None,
    None
]:
    header: Header
    read: AlignedSegment
    results: List[PosNA]
    pbar: Optional[tqdm] = None

    if progress:
        total: int = int(
            pysam.idxstats(samfile)
            .splitlines()[0]
            .split('\t')[2]
        )
        pbar = tqdm(total=total)
    with pysam.AlignmentFile(samfile, 'rb') as samfp:
        for read in samfp.fetch():
            if pbar:
                pbar.update(1)
            if not read.query_sequence:
                continue
            results = iter_single_read_posnas(
                read.query_sequence,
                read.query_qualities,
                read.get_aligned_pairs(False)
            )
            if site_quality_cutoff > 0:
                results = [
                    (pos, idx, na, q)
                    for pos, idx, na, q in results
                    if q >= site_quality_cutoff
                ]
            yield read.query_name, results
        if pbar:
            pbar.close()
