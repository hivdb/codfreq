from tqdm import tqdm  # type: ignore
from array import array
from pysam import AlignedSegment  # type: ignore
from typing import List, Tuple, Optional, Generator, Union

from .codfreq_types import NAPos, NAChar, SeqText, Header
from .paired_reads import PairedReads

ENCODING: str = 'UTF-8'
GAP: int = ord(b'-')

__all__ = ['PosNA', 'iter_single_read_posnas', 'iter_posnas']

PosNA = Tuple[
    Tuple[
        NAPos,  # refpos
        int   # insertion_index
    ],
    NAChar,  # na
    int   # qua
]


def iter_single_read_posnas(
    seq: SeqText,
    qua: array,
    aligned_pairs: List[Tuple[Optional[NAPos], Optional[NAPos]]]
) -> Generator[PosNA, None, None]:
    seqpos0: Optional[NAPos]
    refpos0: Optional[NAPos]
    refpos: NAPos
    n: NAChar  # na
    q: int  # qua
    posna: PosNA

    seqchars: bytearray = bytearray(seq, ENCODING)

    prev_refpos: int = 0
    prev_seqpos0: int = 0
    idx: int = 0
    prev_ins_buffer: List[PosNA] = []

    for seqpos0, refpos0 in aligned_pairs:

        if refpos0 is None:
            # insertion
            refpos = prev_refpos
            idx += 1
        else:
            refpos = refpos0 + 1
            idx = 0
            prev_refpos = refpos

        if seqpos0 is None:
            # deletion
            n, q = GAP, qua[prev_seqpos0]
        else:
            n, q = seqchars[seqpos0], qua[seqpos0]
            prev_seqpos0 = seqpos0

        if refpos == 0:
            continue

        posna = (refpos, idx), n, q
        if idx > 0:
            prev_ins_buffer.append(posna)
        else:
            if prev_ins_buffer:
                yield from prev_ins_buffer
                prev_ins_buffer = []
            yield posna


def iter_posnas(
    all_paired_reads: List[PairedReads],
    site_quality_cutoff: int = 0,
    progress: bool = True
) -> Generator[
    Tuple[str, Generator[PosNA, None, None]],
    None,
    None
]:
    header: Header
    pair: List[PairedReads]
    read: AlignedSegment
    results: Generator[PosNA, None, None]
    _all_paired_reads: Union[List[PairedReads], tqdm]

    if progress:
        _all_paired_reads = tqdm(all_paired_reads)
    else:
        _all_paired_reads = all_paired_reads
    for header, pair in _all_paired_reads:
        for read in pair:
            if not read.query_sequence:
                continue
            results = iter_single_read_posnas(
                read.query_sequence,
                read.query_qualities,
                read.get_aligned_pairs(False)
            )
            if site_quality_cutoff > 0:
                results = (
                    (posidx, na, q)
                    for posidx, na, q in results
                    if q >= site_quality_cutoff
                )
            yield header, results
