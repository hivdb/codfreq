import pysam  # type: ignore
import cython  # type: ignore
from tqdm import tqdm  # type: ignore
from collections import defaultdict, Counter

from typing import (
    Generator,
    Tuple,
    List,
    Dict,
    Optional,
    Any,
    Union,
    DefaultDict,
    Counter as tCounter
)
from concurrent.futures import ProcessPoolExecutor

from .codfreq_types import (
    Header,
    Profile,
    CodFreqRow,
    CodonText,
    GeneInterval,
    FASTQFileName,
    FragmentConfig,
    MainFragmentConfig,
    DerivedFragmentConfig,
)
from .sam2codfreq_types import (
    CodonCounter,
    TypedRefFragment
)
from .samfile_helper import chunked_samfile
from .json_progress import JsonProgress
from .poscodons import iter_poscodons, PosCodon
from .codonalign_consensus import codonalign_consensus
from .filename_helper import name_bamfile


CODFREQ_HEADER: List[str] = [
    'gene', 'position',
    'total', 'codon', 'count',
    'total_quality_score'
]

ENCODING: str = 'UTF-8'


@cython.cfunc
@cython.inline
@cython.returns(list)
def build_gene_intervals(
    genes: List[DerivedFragmentConfig]
) -> List[GeneInterval]:
    return [(
        genedef['refStart'],
        genedef['refEnd'],
        genedef['geneName']
    ) for genedef in genes]


@cython.cfunc
@cython.inline
@cython.returns(list)
def get_ref_fragments(
    profile: Profile
) -> List[
    Tuple[str, MainFragmentConfig, List[DerivedFragmentConfig]]
]:
    refname: str
    refseq: Optional[str]
    fromref: Optional[str]
    refstart: Optional[int]
    refend: Optional[int]
    config: FragmentConfig
    ref_fragments: Dict[str, TypedRefFragment] = {}
    for config in profile['fragmentConfig']:
        refname = config['fragmentName']
        refseq = config.get('refSequence')
        if refseq:
            ref_fragments[refname] = {
                'ref': {
                    'fragmentName': refname,
                    'refSequence': refseq
                },
                'genes': []
            }
    for config in profile['fragmentConfig']:
        refname = config['fragmentName']
        fromref = config.get('fromFragment')
        gene = config.get('geneName')
        refstart = config.get('refStart')
        refend = config.get('refEnd')

        if fromref and gene and refstart and refend:
            ref_fragments[fromref]['genes'].append({
                'fragmentName': refname,
                'fromFragment': fromref,
                'geneName': gene,
                'refStart': refstart,
                'refEnd': refend
            })

    return [
        (refname, pair['ref'], pair['genes'])
        for refname, pair in ref_fragments.items()
    ]


@cython.cfunc
@cython.inline
@cython.returns(list)
def get_codonfreq(
    codonstat: CodonCounter,
    qualities: CodonCounter
) -> List[CodFreqRow]:
    gene: str
    refpos: int
    codons: tCounter[CodonText]
    codon: CodonText
    count: int
    qua: int
    total: int
    rows: List[CodFreqRow] = []
    reformed: DefaultDict[
        Tuple[str, int],
        Counter[CodonText]
    ] = defaultdict(Counter)

    for (gene, refpos, codon), count in codonstat.items():
        reformed[(gene, refpos)][codon] = count

    for (gene, refpos), codons in reformed.items():
        total = sum(codons.values())
        for codon, count in sorted(codons.most_common()):
            qua = qualities[(gene, refpos, codon)]
            rows.append({
                'gene': gene,
                'position': refpos,
                'total': total,
                'codon': codon,
                'count': count,
                'total_quality_score': round(qua, 2)
            })
    return rows


@cython.ccall
@cython.returns(tuple)
def sam2codfreq_between(
    samfile: str,
    samfile_start: int,
    samfile_end: int,
    gene_intervals: List[GeneInterval],
    site_quality_cutoff: int = 0
) -> Tuple[CodonCounter, CodonCounter, int]:
    """subprocess function to call iter_poscodons and count codons

    The results of PosCodons are aggregated into codon counters to reduce
    memory usage before they are sent to the main process. This also reduces
    the computational complexity of the main process and relieves the risk of
    backpressure, which may lead to an out-of-memory issue.

    The out-of-memory issue was misidentified as a mysterious memory leak
    issue. It got worse when Cython was enabled. It turns out that Cython only
    speeded up the subprocesses and sent data to main process in shorter times.
    The main process was simply not fast enough to process such amount of huge
    data. The unprocessed data was piled in the memory of main process and
    caused out-of-memory eventually.
    """

    poscodons: List[PosCodon]
    gene: str
    refpos: int
    codon: CodonText
    qua: int
    codonstat: CodonCounter = Counter()
    qualities: CodonCounter = Counter()
    num_row: int = 0

    for _, poscodons in iter_poscodons(
        samfile,
        samfile_start,
        samfile_end,
        gene_intervals,
        site_quality_cutoff
    ):
        num_row += 1
        for gene, refpos, codon, qua in poscodons:
            codonstat[(gene, refpos, codon)] += 1
            qualities[(gene, refpos, codon)] += qua

    return codonstat, qualities, num_row


def sam2codfreq(
    samfile: str,
    ref: MainFragmentConfig,
    genes: List[DerivedFragmentConfig],
    workers: int,
    site_quality_cutoff: int = 0,
    log_format: str = 'text',
    include_partial_codons: bool = False,
    chunk_size: int = 25000,
    **extras: Any
) -> List[CodFreqRow]:
    """Returns CodFreq rows from a SAM/BAM file

    This function utilizes subprocesses to process alignment data from a
    segment of SAM/BAM file and to aggregate the results into a codon counter
    and a quality score counter.

    The codon counter and quality score counter are converted into CodFreq
    results. The consensus codon from CodFreq results are then processed by
    codonalign function provided by hivdb/postalign.

    :param samfile: str of the SAM/BAM file path
    :param ref: dict of the reference configuration
    :param genes: list of dict of the gene configurations
    :param site_quality_cutoff: phred-score quality cutoff of each codon
                                position
    :param log_format: 'json' or 'text', default to 'text'
    :param include_partial_codons: if true output partial codons also; only for
                                   debug purpose
    :param **extras: any other variables to pass to the log method
    :return: CodFreq rows
    """

    poscodons: Tuple[Header, List[PosCodon]]
    gene: str
    refpos: int
    codon: CodonText
    qua: int
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
    gene_intervals: List[GeneInterval] = build_gene_intervals(genes)
    codonstat: CodonCounter = Counter()
    qualities: CodonCounter = Counter()

    with ProcessPoolExecutor(workers) as executor:

        for partial_codonstat, partial_qualities, num_row in executor.map(
            sam2codfreq_between,
            *zip(*[
                (
                    samfile,
                    samfile_begin,
                    samfile_end,
                    gene_intervals,
                    site_quality_cutoff
                )
                for samfile_begin, samfile_end in chunks
            ])
        ):
            codonstat.update(partial_codonstat)
            qualities.update(partial_qualities)
            if pbar:
                pbar.update(num_row)
        if pbar:
            pbar.close()

    codonfreq: List[CodFreqRow] = get_codonfreq(codonstat, qualities)

    # Apply codon-alignment to consensus codons. This step is an imperfect
    # approach to address the out-frame deletions/insertions. However, it
    # is faster than codon-align each single read, since the process is
    # right now very slow and time-consuming.  This approach may need to be
    # changed in the future when optimization was done for the post-align
    # codon-alignment program.
    codonfreq = codonalign_consensus(codonfreq, ref, genes)
    return codonfreq


def sam2codfreq_all(
    name: str,
    fnpair: Tuple[Optional[FASTQFileName], ...],
    profile: Profile,
    workers: int,
    site_quality_cutoff: int = 0,
    log_format: str = 'text',
    include_partial_codons: bool = False
) -> Generator[CodFreqRow, None, None]:
    refname: str
    ref: MainFragmentConfig
    genes: List[DerivedFragmentConfig]
    for refname, ref, genes in get_ref_fragments(profile):
        samfile: str = name_bamfile(name, refname)
        yield from sam2codfreq(
            samfile,
            ref,
            genes,
            workers,
            site_quality_cutoff=site_quality_cutoff,
            log_format=log_format,
            include_partial_codons=include_partial_codons,
            fastqs=fnpair
        )
