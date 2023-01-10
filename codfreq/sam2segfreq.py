import pysam  # type: ignore
import cython  # type: ignore
from tqdm import tqdm  # type: ignore

from typing import (
    Tuple,
    List,
    Optional,
    Any,
    Union,
    Deque
)
from concurrent.futures import ProcessPoolExecutor
from collections import deque

from .codfreq_types import (
    Profile,
    FASTQFileName,
    MainFragmentConfig
)
from .samfile_helper import chunked_samfile
from .json_progress import JsonProgress
from .posnas import PosNA, iter_single_read_posnas
from .filename_helper import name_bamfile, name_segfreq
from .profile import get_ref_fragments

from .segfreq import SegFreq, DEFAULT_SEGMENT_SIZE, DEFAULT_SEGMENT_STEP

BEGIN: Tuple[int, int, int, int] = (0, 0, ord('^'), 0)
END: Tuple[int, int, int, int] = (0, 0, ord('$'), 0)
ENCODING: str = 'UTF-8-sig'


@cython.ccall
@cython.locals(
    samfile=str,
    samfile_start=cython.long,
    samfile_end=cython.long,
    segment_size=cython.int,
    segment_step=cython.int
)
@cython.returns(tuple)
def sam2segfreq_between(
    samfile: str,
    samfile_start: int,
    samfile_end: int,
    segment_size: int,
    segment_step: int
) -> Tuple[SegFreq, int]:
    """subprocess function to count PosNA and edges

    The results of PosNAs are aggregated into edge counters to reduce memory
    usage before they are sent back to the main process.
    """
    fragment_name: str
    segment_deque: Deque[Tuple[int, Optional[PosNA]]]
    segment_pos: int = cython.declare(cython.long)
    next_segment_pos: int = cython.declare(cython.long)
    segfreq: SegFreq = SegFreq(segment_size, segment_step)
    idx: int = cython.declare(cython.int)
    pos: int = cython.declare(cython.long)
    last_pos: int = cython.declare(cython.long)
    num_row: int = cython.declare(cython.long, 0)

    with pysam.AlignmentFile(samfile, 'rb') as samfp:
        samfp.seek(samfile_start)

        for read in samfp:
            num_row += 1
            if samfp.tell() > samfile_end:
                break

            if not read.query_sequence:
                continue

            posnas: List[Tuple[int, Optional[PosNA]]] = [
                (posna.pos, posna)
                for posna in iter_single_read_posnas(
                    read.query_sequence,
                    read.get_aligned_pairs(matches_only=False)
                )
            ]
            # add trailing None's to indicate missing of positions
            last_pos = posnas[-1][0]
            post_posnas: List[Tuple[int, Optional[PosNA]]] = [
                (pos, None)
                for pos in range(last_pos + 1, last_pos + segment_size - 1)
            ]
            posnas = posnas + post_posnas

            segment_pos = 0
            segment_deque = deque()
            for idx, (pos, posna) in enumerate(posnas):

                if segment_pos > 0:
                    size = pos - segment_pos
                    if size > segment_size:
                        raise ValueError(
                            'The read alignment is not continuous')
                    elif size == segment_size:
                        segfreq.add(
                            tuple([posna for _, posna in segment_deque]))
                        next_segment_pos = segment_pos + segment_step
                        while segment_deque[0][0] < next_segment_pos:
                            segment_deque.popleft()
                        segment_pos = next_segment_pos
                elif segment_pos == 0:
                    remainder = (pos - 1) % segment_step
                    if remainder == 0:
                        segment_pos = pos
                    else:
                        segment_pos = pos - remainder
                        segment_deque.extend([
                            (p, None) for p in range(segment_pos, pos)
                        ])

                segment_deque.append((pos, posna))
            else:
                if segment_pos > 0:
                    size = pos - segment_pos + 1
                    if size > segment_size:
                        raise ValueError(
                            'The read alignment is not continuous')
                    elif size == segment_size:
                        segfreq.add(
                            tuple([posna for _, posna in segment_deque]))
    return segfreq, num_row


def sam2segfreq(
    samfile: str,
    ref: MainFragmentConfig,
    workers: int,
    segment_size: int = DEFAULT_SEGMENT_SIZE,
    segment_step: int = DEFAULT_SEGMENT_STEP,
    log_format: str = 'text',
    chunk_size: int = 25000,
    **extras: Any
) -> SegFreq:
    """Returns a SegFreq object from a SAM/BAM file

    This function utilizes subprocesses to process alignment data from a
    segment of SAM/BAM file and to aggregate the results into an DiGraph edge
    counter (SegFreq).

    :param samfile: str of the SAM/BAM file path
    :param ref: dict of the reference configuration
    :param workers: int of the number of subprocesses
    :param segment_size: int of the segment size, default 3
    :param segment_step: int of the segment step, default 1
    :param log_format: 'json' or 'text', default to 'text'
    :param **extras: any other variables to pass to the log method
    :return: a SegFreq object
    """

    total: int
    pbar: Optional[Union[JsonProgress, tqdm]]

    with pysam.AlignmentFile(samfile, 'rb') as samfp:
        total = samfp.mapped
        if log_format == 'json':
            pbar = JsonProgress(
                total=total, description=samfile, **extras)
        elif log_format == 'text':
            pbar = tqdm(total=total)
            pbar.set_description('Processing {}'.format(samfile))

    chunks: List[Tuple[int, int]] = chunked_samfile(samfile, chunk_size)
    segfreq: SegFreq = SegFreq(segment_size, segment_step)

    with ProcessPoolExecutor(workers) as executor:

        for partial_segfreq, num_row in executor.map(
            sam2segfreq_between,
            *zip(*[
                (
                    samfile,
                    samfile_begin,
                    samfile_end,
                    segment_size,
                    segment_step
                )
                for samfile_begin, samfile_end in chunks
            ])
        ):
            segfreq.update(partial_segfreq)
            if pbar:
                pbar.update(num_row)
        if pbar:
            pbar.close()

    return segfreq


def sam2segfreq_all(
    name: str,
    fnpair: Tuple[Optional[FASTQFileName], ...],
    profile: Profile,
    workers: int,
    log_format: str = 'text'
) -> None:
    refname: str
    ref: MainFragmentConfig
    ref_fragments = get_ref_fragments(profile)
    for refname, ref, _ in ref_fragments:
        samfile: str = name_bamfile(name, refname, is_trimmed=True)
        segfreq: SegFreq = sam2segfreq(
            samfile,
            ref,
            workers=workers,
            segment_size=ref['segmentSize'],
            segment_step=ref['segmentStep'],
            log_format=log_format,
            fastqs=fnpair
        )
        with open(
            name_segfreq(name, refname),
            'w',
            encoding=ENCODING
        ) as segfreqfp:
            segfreq.dump(segfreqfp)
