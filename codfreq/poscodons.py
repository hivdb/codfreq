import pysam  # type: ignore
from tqdm import tqdm  # type: ignore
from itertools import groupby
from operator import itemgetter
from pysam import AlignedSegment  # type: ignore
from typing import List, Tuple, Any, Union, Generator, Iterator

from .codfreq_types import (
    MainFragmentConfig,
    DerivedFragmentConfig,
    Header,
    AAPos,
    NAPos,
    GeneText,
    CodonText
)
from .json_progress import JsonProgress
from .posnas import iter_single_read_posnas, PosNA

# frameshift types
REPEAT: int = 0x01
SKIP: int = 0x02

UNSEQ: bytes = b'.'

#                  refStart refEnd
#                      v      v
GeneInterval = Tuple[NAPos, NAPos, GeneText]
#                                                        Qual
#                                                         v
PosCodon = Tuple[GeneText, AAPos, Tuple[CodonText, ...], int]
GroupedPosNAs = List[
    Tuple[GeneText, AAPos, List[PosNA]]
]


def build_gene_intervals(
    genes: List[DerivedFragmentConfig]
) -> List[GeneInterval]:
    return [(
        genedef['refStart'],
        genedef['refEnd'],
        genedef['geneName']
    ) for genedef in genes]


def pos_from_posna(posna: PosNA) -> int:
    return posna[0][0]


def group_posnas(
    posnas: List[PosNA],
    gene_intervals: List[GeneInterval]
) -> GroupedPosNAs:
    """Group posnas by position and add gene AA position"""

    posnas_by_napos: List[Tuple[NAPos, List[PosNA]]] = [
        (pos, list(na_and_ins))
        for pos, na_and_ins in
        groupby(posnas, key=pos_from_posna)
    ]

    pos: NAPos
    na_and_ins: List[PosNA]
    aapos: AAPos
    gene_refstart: NAPos
    gene_refend: NAPos
    gene: GeneText
    curidx: int = -1
    size: int = len(posnas_by_napos)
    grouped_posnas: GroupedPosNAs = []

    # convert napos to gene and aapos
    for gene_refstart, gene_refend, gene in gene_intervals:
        pos = 0
        # seek back if pos at curidx is after gene_refstart.
        # this deals with situations such as the overlap of
        # ORF1a and ORF1b.
        while (
            curidx > -1 and
            curidx + 1 < size and
            posnas_by_napos[curidx + 1][0] > gene_refstart
        ):
            curidx -= posnas_by_napos[curidx + 1][0] - gene_refstart
        # The above while loop can set curidx < -1. We had to change
        # it back to -1 or the index applied to posnas_by_napos below
        # will be a negative number.
        if curidx < -1:
            curidx = -1
        while curidx + 1 < size and pos < gene_refend:
            curidx += 1
            # This number should be zero or positive
            #                                   v
            pos, na_and_ins = posnas_by_napos[curidx]
            if pos < gene_refstart:
                continue
            aapos = (pos - gene_refstart) // 3 + 1

            grouped_posnas.append((
                gene,
                aapos,
                na_and_ins
            ))
    return grouped_posnas


def find_overlapped_genes(
    gene_intervals: List[GeneInterval],
    read_refstart: NAPos,
    read_refend: NAPos
) -> List[GeneInterval]:
    gene_refstart: NAPos
    gene_refend: NAPos
    gene: GeneText
    filtered: List[GeneInterval] = []
    for gene_refstart, gene_refend, gene in gene_intervals:
        if read_refend < gene_refstart:
            break
        if read_refstart > gene_refend:
            continue
        filtered.append((gene_refstart, gene_refend, gene))
    return filtered


def get_comparable_codon(
    codon_posnas: List[List[PosNA]]
) -> Tuple[Tuple[CodonText, ...], bool]:
    posnas: List[PosNA]
    codon_nas: Tuple[CodonText, ...] = tuple([
        bytes([na for _, na, _ in posnas])
        for posnas in codon_posnas
    ])
    is_partial: bool = len(codon_nas) < 3
    return codon_nas, is_partial


def posnas2poscodons(
    posnas: List[PosNA],
    gene_intervals: List[GeneInterval],
    read_refstart: int,  # 1-based first aligned refpos
    read_refend: int,    # 1-based last aligned refpos
    site_quality_cutoff: int
) -> List[PosCodon]:
    meanq: List[int]
    meanq_int: int
    gene: GeneText
    aapos: AAPos
    codon_iter: Iterator[Tuple[GeneText, AAPos, List[PosNA]]]
    codon_posnas: List[List[PosNA]]
    codons: Tuple[CodonText, ...]
    is_partial: bool

    genes: List[GeneInterval] = find_overlapped_genes(
        gene_intervals, read_refstart, read_refend)
    grouped_posnas: GroupedPosNAs = group_posnas(posnas, genes)

    poscodons: List[PosCodon] = []
    for (gene, aapos), codon_iter in groupby(
        grouped_posnas, key=itemgetter(0, 1)
    ):
        codon_posnas = [pnas for _, _, pnas in codon_iter]
        meanq = [qua for pnas in codon_posnas for _, _, qua in pnas]
        meanq_int = round(sum(meanq) / len(meanq) if meanq else 0)
        if meanq_int < site_quality_cutoff:
            continue
        codons, is_partial = get_comparable_codon(codon_posnas)
        if is_partial:
            continue
        poscodons.append((gene, aapos, codons, meanq_int))
    return poscodons


def iter_poscodons(
    samfile: str,
    ref: MainFragmentConfig,
    genes: List[DerivedFragmentConfig],
    description: str,
    site_quality_cutoff: int = 0,
    log_format: str = 'text',
    **extras: Any
) -> Generator[
    Tuple[Header, List[PosCodon]],
    None,
    None
]:

    pbar: Union[JsonProgress, tqdm]
    header: Header
    read: AlignedSegment
    pair: List[AlignedSegment]
    posnas: List[PosNA]
    poscodons: List[PosCodon]
    gene_intervals: List[GeneInterval] = build_gene_intervals(genes)
    total: int = int(
        pysam.idxstats(samfile)
        .splitlines()[0]
        .split('\t')[2]
    )
    if log_format == 'json':
        pbar = JsonProgress(
            total=total, description=description, **extras)
    else:
        pbar = tqdm(total=total)
        pbar.set_description('Processing {}'.format(description))

    with pysam.AlignmentFile(samfile, 'rb') as samfp:
        for read in samfp.fetch():
            pbar.update(1)
            if not read.query_sequence:
                continue

            posnas = iter_single_read_posnas(
                read.query_sequence,
                read.query_qualities,
                read.get_aligned_pairs(False)
            )

            poscodons = posnas2poscodons(
                posnas,
                gene_intervals,
                read.reference_start + 1,  # pysam has 0-based numbering
                read.reference_end,  # "reference_end points to one past the
                                     #  last aligned residue."
                site_quality_cutoff
            )
            yield read.query_name, poscodons
        pbar.close()
