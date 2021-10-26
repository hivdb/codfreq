from tqdm import tqdm  # type: ignore
from itertools import groupby
from operator import itemgetter
from pysam import AlignedSegment  # type: ignore
from typing import List, Tuple, Any, Union, Generator, Iterator

from .paired_reads import PairedReads
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

RefStart = NAPos
RefEnd = NAPos
Qua = int
GeneInterval = Tuple[RefStart, RefEnd, GeneText]
PosCodon = Tuple[GeneText, AAPos, Tuple[CodonText, ...], Qua]
GroupedPosNAs = Generator[
    Tuple[GeneText, AAPos, List[PosNA]], None, None
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
    posnas: Generator[PosNA, None, None],
    gene_intervals: Generator[GeneInterval, None, None]
) -> GroupedPosNAs:
    """Group posnas by position and add gene AA position"""

    grouped_posnas: List[Tuple[int, List[PosNA]]] = [
        (pos, list(na_and_ins))
        for pos, na_and_ins in
        groupby(posnas, key=pos_from_posna)
    ]

    pos: NAPos
    na_and_ins: List[PosNA]
    gene_refstart: RefStart
    gene_refend: RefEnd
    gene: GeneText
    curidx: int = -1
    size: int = len(grouped_posnas)
    for gene_refstart, gene_refend, gene in gene_intervals:
        pos = 0
        # seek back if pos at curidx is after gene_refstart
        while (
            curidx > -1 and
            curidx + 1 < size and
            grouped_posnas[curidx + 1][0] > gene_refstart
        ):
            curidx -= grouped_posnas[curidx + 1][0] - gene_refstart
        while curidx + 1 < size and pos < gene_refend:
            curidx += 1
            pos, na_and_ins = grouped_posnas[curidx]
            if pos < gene_refstart:
                continue
            aapos = (pos - gene_refstart) // 3 + 1
            # if gene == 'nsp1' and aapos == 90:
            #     print(gene, aapos, na_and_ins_tuple)

            yield (
                gene,
                aapos,
                na_and_ins
            )


def find_overlapped_genes(
    gene_intervals: List[GeneInterval],
    read_refstart: RefStart,
    read_refend: RefEnd
) -> Generator[GeneInterval, None, None]:
    gene_refstart: RefStart
    gene_refend: RefEnd
    gene: GeneText
    for gene_refstart, gene_refend, gene in gene_intervals:
        if read_refend < gene_refstart:
            break
        if read_refstart > gene_refend:
            continue
        yield gene_refstart, gene_refend, gene


def get_comparable_codon(
    codon_posnas: List[List[PosNA]]
) -> Tuple[Tuple[CodonText, ...], bool]:
    posnas: List[PosNA]
    codon_nas: Tuple[CodonText, ...] = tuple(
        bytes(na for _, na, _ in posnas)
        for posnas in codon_posnas
    )
    is_partial: bool = len(codon_nas) < 3
    return codon_nas, is_partial


def posnas2poscodons(
    posnas: Generator[PosNA, None, None],
    gene_intervals: List[GeneInterval],
    read_refstart: int,  # 1-based first aligned refpos
    read_refend: int,    # 1-based last aligned refpos
    site_quality_cutoff: int
) -> Generator[PosCodon, None, None]:
    meanq: List[Qua]
    meanq_int: Qua
    gene: GeneText
    aapos: AAPos
    codon_iter: Iterator[Tuple[GeneText, AAPos, List[PosNA]]]
    codon_posnas: List[List[PosNA]]
    codons: Tuple[CodonText, ...]
    is_partial: bool

    genes: Generator[GeneInterval, None, None] = find_overlapped_genes(
        gene_intervals, read_refstart, read_refend)
    grouped_posnas: GroupedPosNAs = group_posnas(posnas, genes)
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
        yield gene, aapos, codons, meanq_int


def iter_poscodons(
    all_paired_reads: List[PairedReads],
    ref: MainFragmentConfig,
    genes: List[DerivedFragmentConfig],
    description: str,
    site_quality_cutoff: int = 0,
    log_format: str = 'text',
    **extras: Any
) -> Generator[
    Tuple[Header, Generator[PosCodon, None, None]],
    None,
    None
]:
    pbar: Union[JsonProgress, tqdm]
    header: Header
    read: AlignedSegment
    pair: List[AlignedSegment]
    posnas: Generator[PosNA, None, None]
    poscodons: Generator[PosCodon, None, None]
    gene_intervals: List[GeneInterval] = build_gene_intervals(genes)
    total: int = len(all_paired_reads)

    if log_format == 'json':
        pbar = JsonProgress(
            total=total, description=description, **extras)
    else:
        pbar = tqdm(total=total)
        pbar.set_description('Processing {}'.format(description))

    for header, pair in all_paired_reads:
        pbar.update(1)
        for read in pair:
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
            yield header, poscodons
    pbar.close()
