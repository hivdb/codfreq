import pysam  # type: ignore
import cython  # type: ignore
from tqdm import tqdm  # type: ignore

from typing import (
    Tuple,
    List,
    Dict,
    Optional,
    Any,
    Union
)
from concurrent.futures import ProcessPoolExecutor

from .codfreq_types import (
    Header,
    Profile,
    FASTQFileName,
    FragmentConfig,
    MainFragmentConfig
)
from .sam2codfreq_types import (
    TypedRefFragment
)
from .samfile_helper import chunked_samfile
from .json_progress import JsonProgress
from .posnas import PosNA, iter_single_read_posnas
from .filename_helper import name_bamfile, name_segfreq

from .segfreq import SegFreq

BEGIN: Tuple[int, int, int, int] = (0, 0, ord('^'), 0)
END: Tuple[int, int, int, int] = (0, 0, ord('$'), 0)
ENCODING: str = 'UTF-8-sig'


@cython.cfunc
@cython.inline
@cython.returns(list)
def get_ref_fragments(
    profile: Profile
) -> List[Tuple[
    Header,
    MainFragmentConfig
]]:
    refname: str
    config: FragmentConfig
    ref_fragments: Dict[Header, TypedRefFragment] = {}
    for config in profile['fragmentConfig']:
        refname = config['fragmentName']
        refseq = config.get('refSequence')
        if isinstance(refseq, str):
            ref_fragments[refname] = {
                'ref': {
                    'fragmentName': refname,
                    'refSequence': refseq
                },
                'fragments': []
            }

    return [
        (refname, pair['ref'])
        for refname, pair in ref_fragments.items()
    ]


@cython.ccall
@cython.locals(
    samfile=str,
    samfile_start=cython.ulong,
    samfile_end=cython.ulong,
    segment_size=cython.uint,
    segment_step=cython.uint
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
    segfreq: SegFreq = SegFreq(segment_size, segment_step)
    idx: int = cython.declare(cython.uint)
    jdx: int = cython.declare(cython.uint)
    pos: int = cython.declare(cython.ulong)
    first_pos: int = cython.declare(cython.ulong)
    last_pos: int = cython.declare(cython.ulong)
    num_row: int = cython.declare(cython.ulong, 0)
    segment_list: List[Optional[PosNA]] = cython.declare(list)

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
            first_pos = posnas[0][0]
            last_pos = posnas[-1][0]
            pre_posnas: List[Tuple[int, Optional[PosNA]]] = [
                (pos, None)
                for pos in range(first_pos - segment_step + 1, first_pos)
            ]
            post_posnas: List[Tuple[int, Optional[PosNA]]] = [
                (pos, None)
                for pos in range(last_pos + 1, last_pos + segment_step - 1)
            ]
            posnas = pre_posnas + posnas + post_posnas

            for idx, (first_pos, first_posna) in enumerate(posnas):
                if (first_pos - 1) % segment_step != 0:
                    break
                elif first_posna and first_posna.bp > 0:
                    break
                segment_list = [first_posna]

                posnas_len = len(posnas)
                jdx = idx + 1
                while jdx < posnas_len:
                    pos, next_posna = posnas[jdx]
                    jdx += 1
                    segment_list.append(next_posna)
                    if pos - first_pos + 1 == segment_size:
                        segfreq.add(tuple(segment_list))
                        break
    return segfreq, num_row


def sam2segfreq(
    samfile: str,
    ref: MainFragmentConfig,
    workers: int,
    segment_size: int = 12,
    segment_step: int = 3,
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
    :param segment_size: int of the segment size, default 12
    :param segment_step: int of the segment step, default 3
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
    log_format: str = 'text',
    include_partial_codons: bool = False
) -> None:
    refname: str
    ref: MainFragmentConfig
    ref_fragments = get_ref_fragments(profile)
    for refname, ref in ref_fragments:
        samfile: str = name_bamfile(name, refname, is_trimmed=True)
        segfreq: SegFreq = sam2segfreq(
            samfile,
            ref,
            workers=workers,
            segment_size=profile.get('segmentSize', 6),
            segment_step=profile.get('segmentStep', 1),
            log_format=log_format,
            fastqs=fnpair
        )
        with open(
            name_segfreq(name, refname),
            'w',
            encoding=ENCODING
        ) as segfreqfp:
            segfreq.dump(segfreqfp)
