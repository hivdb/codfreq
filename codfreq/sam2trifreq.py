import pysam  # type: ignore
import cython  # type: ignore
from tqdm import tqdm  # type: ignore
from more_itertools import triplewise

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
from .posnas import iter_single_read_posnas
from .filename_helper import name_bamfile, name_trifreq

from .trifreq import PosNA, TriFreq

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
@cython.returns(tuple)
def sam2trifreq_between(
    samfile: str,
    samfile_start: int,
    samfile_end: int
) -> Tuple[TriFreq, int]:
    """subprocess function to count PosNA and edges

    The results of PosNAs are aggregated into edge counters to reduce memory
    usage before they are sent back to the main process.
    """
    fragment_name: str
    trifreq: TriFreq = TriFreq()
    num_row: int = 0

    with pysam.AlignmentFile(samfile, 'rb') as samfp:
        samfp.seek(samfile_start)

        for read in samfp:
            num_row += 1
            if samfp.tell() > samfile_end:
                break

            if not read.query_sequence:
                continue

            posnas = iter_single_read_posnas(
                read.query_sequence,
                read.query_qualities,
                read.get_aligned_pairs(matches_only=False)
            )

            for (
                pos_a, bp_a, na_a, _
            ), (
                pos_b, bp_b, na_b, _
            ), (
                pos_c, bp_c, na_c, _
            ) in triplewise(posnas):
                posna_a = PosNA(pos_a, bp_a, na_a)
                posna_b = PosNA(pos_b, bp_b, na_b)
                posna_c = PosNA(pos_c, bp_c, na_c)
                trifreq.add_triplet(posna_a, posna_b, posna_c)
    return trifreq, num_row


def sam2trifreq(
    samfile: str,
    ref: MainFragmentConfig,
    workers: int,
    log_format: str = 'text',
    chunk_size: int = 25000,
    **extras: Any
) -> TriFreq:
    """Returns a TriFreq object from a SAM/BAM file

    This function utilizes subprocesses to process alignment data from a
    segment of SAM/BAM file and to aggregate the results into an DiGraph edge
    counter (TriFreq).

    :param samfile: str of the SAM/BAM file path
    :param ref: dict of the reference configuration
    :param workers: int of the number of subprocesses
    :param log_format: 'json' or 'text', default to 'text'
    :param **extras: any other variables to pass to the log method
    :return: a TriFreq object
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
    trifreq: TriFreq = TriFreq()

    with ProcessPoolExecutor(workers) as executor:

        for partial_trifreq, num_row in executor.map(
            sam2trifreq_between,
            *zip(*[
                (
                    samfile,
                    samfile_begin,
                    samfile_end
                )
                for samfile_begin, samfile_end in chunks
            ])
        ):
            trifreq.update(partial_trifreq)
            if pbar:
                pbar.update(num_row)
        if pbar:
            pbar.close()

    return trifreq


def sam2trifreq_all(
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
        trifreq: TriFreq = sam2trifreq(
            samfile,
            ref,
            workers=workers,
            log_format=log_format,
            fastqs=fnpair
        )
        with open(
            name_trifreq(name, refname),
            'w',
            encoding=ENCODING
        ) as trifreqfp:
            trifreq.dump(trifreqfp)
