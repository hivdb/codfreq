import cython  # type: ignore
from collections import defaultdict, Counter

import typing
from typing import (
    Generator,
    Tuple,
    List,
    Dict,
    Optional,
    Any
)
from .codfreq_types import (
    Profile,
    CodFreqRow,
    MultiNAText,
    FASTQFileName,
    FragmentConfig,
    MainFragmentConfig,
    DerivedFragmentConfig,
)
from .sam2codfreq_types import (
    CodonCounter,
    QualityCounter,
    TypedRefFragment
)
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
    qualities: QualityCounter
) -> List[CodFreqRow]:
    gene: str
    refpos: int
    codons: typing.Counter[Tuple[MultiNAText, ...]]
    codon: Tuple[MultiNAText, ...]
    count: int
    total: int
    qua: int
    rows: List[CodFreqRow] = []

    for (gene, refpos), codons in codonstat.items():
        total = sum(codons.values())
        for codon, count in sorted(codons.most_common()):
            qua = qualities[(gene, refpos, codon)]
            rows.append({
                'gene': gene,
                'position': refpos,
                'total': total,
                'codon': bytes([na for bp in codon for na in bp]),
                'count': count,
                'total_quality_score': round(qua, 2)
            })
    return rows


def sam2codfreq(
    samfile: str,
    ref: MainFragmentConfig,
    genes: List[DerivedFragmentConfig],
    workers: int,
    site_quality_cutoff: int = 0,
    log_format: str = 'text',
    include_partial_codons: bool = False,
    **extras: Any
) -> List[CodFreqRow]:
    """Iterate CodFreq rows from a SAM/BAM file

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

    poscodons: List[PosCodon]
    gene: str
    refpos: int
    codon: Tuple[bytes, ...]
    qua: int
    codonstat: CodonCounter = defaultdict(Counter)
    qualities: QualityCounter = Counter()
    # Iterate through the whole SAM/BAM and stat each individual gene
    # position and codon
    for _, poscodons in iter_poscodons(
        samfile,
        genes,
        workers,
        description=samfile,
        site_quality_cutoff=site_quality_cutoff,
        log_format=log_format,
        include_partial_codons=include_partial_codons,
        **extras
    ):
        for gene, refpos, codon, qua in poscodons:
            codonstat[(gene, refpos)][codon] += 1
            qualities[(gene, refpos, codon)] += qua
        del poscodons

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
