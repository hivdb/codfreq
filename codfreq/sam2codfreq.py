from more_itertools import flatten  # type: ignore
from collections import defaultdict, Counter

import typing
from typing import (
    Generator,
    DefaultDict,
    Tuple,
    List,
    Dict,
    TypedDict,
    Optional,
    Any
)
from .codfreq_types import (
    Profile,
    CodFreqRow,
    FASTQFileName,
    FragmentConfig,
    MainFragmentConfig,
    DerivedFragmentConfig,
)
from .paired_reads import iter_paired_reads, PairedReads
from .poscodons import iter_poscodons, PosCodon
from .codonalign_consensus import codonalign_consensus
from .filename_helper import name_bamfile

CODFREQ_HEADER: List[str] = [
    'gene', 'position',
    'total', 'codon', 'count',
    'total_quality_score'
]

ENCODING: str = 'UTF-8'


class TypedRefFragment(TypedDict):
    ref: MainFragmentConfig
    genes: List[DerivedFragmentConfig]


CodonStat = DefaultDict[
    Tuple[str, int],
    typing.Counter[Tuple[bytes, ...]]
]

Qualities = typing.Counter[
    Tuple[str, int, Tuple[bytes, ...]]
]


def iter_ref_fragments(
    profile: Profile
) -> Generator[
    Tuple[str, MainFragmentConfig, List[DerivedFragmentConfig]],
    None,
    None
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

    for refname, pair in ref_fragments.items():
        yield refname, pair['ref'], pair['genes']


def iter_codonfreq(
    codonstat: CodonStat,
    qualities: Qualities
) -> Generator[CodFreqRow, None, None]:
    gene: str
    refpos: int
    codons: typing.Counter[Tuple[bytes, ...]]
    codon: Tuple[bytes, ...]
    count: int
    total: int
    qua: int
    for (gene, refpos), codons in codonstat.items():
        total = sum(codons.values())
        for codon, count in sorted(codons.most_common()):
            qua = qualities[(gene, refpos, codon)]
            yield {
                'gene': gene,
                'position': refpos,
                'total': total,
                'codon': bytes(flatten(codon)),
                'count': count,
                'total_quality_score': round(qua, 2)
            }


def sam2codfreq(
    samfile: str,
    ref: MainFragmentConfig,
    genes: List[DerivedFragmentConfig],
    site_quality_cutoff: int = 0,
    log_format: str = 'text',
    include_partial_codons: bool = False,
    **extras: Any
) -> Generator[CodFreqRow, None, None]:
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

    poscodons: Generator[PosCodon, None, None]
    gene: str
    refpos: int
    codon: Tuple[bytes, ...]
    qua: int
    codonstat: CodonStat = defaultdict(Counter)
    qualities: Qualities = Counter()
    all_paired_reads: List[PairedReads] = iter_paired_reads(samfile)

    # Iterate through the whole SAM/BAM and stat each individual gene position
    # and codon
    for _, poscodons in iter_poscodons(
        all_paired_reads,
        ref,
        genes,
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

    codonfreq: Generator[
        CodFreqRow, None, None
    ] = iter_codonfreq(codonstat, qualities)

    # Apply codon-alignment to consensus codons. This step is an imperfect
    # approach to address the out-frame deletions/insertions. However, it is
    # faster than codon-align each single read, since the process is right
    # now very slow and time-consuming.
    # This approach may need to be changed in the future when optimization was
    # done for the post-align codon-alignment program.
    yield from codonalign_consensus(codonfreq, ref, genes)


def sam2codfreq_all(
    name: str,
    fnpair: Tuple[Optional[FASTQFileName], ...],
    profile: Profile,
    site_quality_cutoff: int = 0,
    log_format: str = 'text',
    include_partial_codons: bool = False
) -> Generator[CodFreqRow, None, None]:
    refname: str
    ref: MainFragmentConfig
    genes: List[DerivedFragmentConfig]
    for refname, ref, genes in iter_ref_fragments(profile):
        samfile: str = name_bamfile(name, refname)
        yield from sam2codfreq(
            samfile,
            ref,
            genes,
            site_quality_cutoff=site_quality_cutoff,
            log_format=log_format,
            include_partial_codons=include_partial_codons,
            fastqs=fnpair
        )
