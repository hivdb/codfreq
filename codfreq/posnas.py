import pysam  # type: ignore
import cython  # type: ignore
from tqdm import tqdm  # type: ignore
from array import array
from pysam import AlignedSegment  # type: ignore
from typing import List, Tuple, Optional, Generator, Any, Union
from concurrent.futures import ProcessPoolExecutor

from .json_progress import JsonProgress
from .samfile_helper import chunked_samfile
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
    qua: Optional[array],
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
            n, q = GAP, qua[prev_seqpos0] if qua else 1
        else:
            n, q = seqchars[seqpos0], qua[seqpos0] if qua else 1
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


@cython.ccall
@cython.returns(list)
def get_posnas_between(
    samfile: str,
    samfile_start: int,
    samfile_end: int,
    site_quality_cutoff: int = 0
) -> List[Tuple[Header, List[PosNA]]]:

    pos: NAPos
    idx: int
    na: NAChar
    q: int
    read: AlignedSegment
    posnas: List[PosNA]

    results: List[Tuple[Header, List[PosNA]]] = []

    with pysam.AlignmentFile(samfile, 'rb') as samfp:

        samfp.seek(samfile_start)

        for read in samfp:
            if samfp.tell() > samfile_end:
                break

            if not read.query_sequence:
                continue

            posnas = iter_single_read_posnas(
                read.query_sequence,
                read.query_qualities,
                read.get_aligned_pairs(False)
            )
            if site_quality_cutoff > 0:
                posnas = [
                    (pos, idx, na, q)
                    for pos, idx, na, q in posnas
                    if q >= site_quality_cutoff
                ]
            results.append((read.query_name, posnas))

    return results


@cython.ccall
@cython.returns(list)
def get_posnas_in_genome_region(
    samfile: str,
    ref_name: str,
    ref_start: NAPos,
    ref_end: NAPos,
    site_quality_cutoff: int = 0
) -> List[Tuple[Header, List[PosNA]]]:

    pos: NAPos
    idx: int
    na: NAChar
    q: int
    read: AlignedSegment
    posnas: List[PosNA]

    results: List[Tuple[Header, List[PosNA]]] = []

    with pysam.AlignmentFile(samfile, 'rb') as samfp:

        for read in samfp.fetch(ref_name, ref_start - 1, ref_end):

            if not read.query_sequence:
                continue

            posnas = iter_single_read_posnas(
                read.query_sequence,
                read.query_qualities,
                read.get_aligned_pairs(False)
            )
            if site_quality_cutoff > 0:
                posnas = [
                    (pos, idx, na, q)
                    for pos, idx, na, q in posnas
                    if q >= site_quality_cutoff
                ]
            results.append((read.query_name, posnas))

    return results


def iter_posnas(
    samfile: str,
    workers: int,
    description: str,
    jsonop: str = 'progress',
    site_quality_cutoff: int = 0,
    chunk_size: int = 5000,
    log_format: str = 'text',
    **extras: Any
) -> Generator[
    Tuple[Header, List[PosNA]],
    None,
    None
]:
    total: int
    header: Header
    chunks: List[Tuple[int, int]]
    samfile_begin: int
    samfile_end: int
    one: Tuple[Header, List[PosNA]]
    posnas: List[Tuple[Header, List[PosNA]]]
    pbar: Optional[Union[JsonProgress, tqdm]] = None

    with pysam.AlignmentFile(samfile, 'rb') as samfp:
        total = samfp.mapped
        if log_format == 'json':
            pbar = JsonProgress(
                op=jsonop,
                total=total,
                description=description,
                **extras)
        elif log_format == 'text':
            pbar = tqdm(total=total)
            pbar.set_description('Processing {}'.format(description))

    chunks = chunked_samfile(samfile, chunk_size)

    with ProcessPoolExecutor(workers) as executor:

        for posnas in executor.map(
            get_posnas_between,
            *zip(*[
                (
                    samfile,
                    samfile_begin,
                    samfile_end,
                    site_quality_cutoff
                )
                for samfile_begin, samfile_end in chunks
            ])
        ):
            if pbar:
                for one in posnas:
                    yield one
                    pbar.update(1)
            else:
                yield from posnas
    if pbar:
        pbar.close()
