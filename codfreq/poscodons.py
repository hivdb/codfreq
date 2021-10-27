import pysam  # type: ignore
import cython  # type: ignore
from tqdm import tqdm  # type: ignore
from pysam import AlignedSegment  # type: ignore
from typing import List, Tuple, Any, Union, Generator, Iterator, Optional
from concurrent.futures import ProcessPoolExecutor

from .codfreq_types import (
    DerivedFragmentConfig,
    Header,
    AAPos,
    NAPos,
    GeneText,
    MultiNAText
)
from .samfile_helper import chunked_samfile
from .json_progress import JsonProgress
from .posnas import iter_single_read_posnas, PosNA

# frameshift types
REPEAT: int = 0x01
SKIP: int = 0x02

UNSEQ: bytes = b'.'

#                  refStart refEnd
#                      v      v
GeneInterval = Tuple[NAPos, NAPos, GeneText]
#                                                          Qual
#                                                           v
PosCodon = Tuple[GeneText, AAPos, Tuple[MultiNAText, ...], int]
BasePair = Tuple[GeneText, AAPos, List[PosNA]]


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
def group_posnas_by_napos(
    posnas: List[PosNA]
) -> List[Tuple[NAPos, List[PosNA]]]:
    prev_pos: int = -1
    by_napos: List[Tuple[NAPos, List[PosNA]]] = []
    for posna in posnas:
        if prev_pos == posna[0]:
            by_napos[-1][1].append(posna)
        else:
            prev_pos = posna[0]
            by_napos.append((posna[0], [posna]))
    return by_napos


@cython.cfunc
@cython.inline
@cython.returns(list)
def group_basepairs(
    posnas: List[PosNA],
    gene_intervals: List[GeneInterval]
) -> List[BasePair]:
    """Group same base-pair posnas by its gene AA position"""

    pos: NAPos
    na_and_ins: List[PosNA]
    aapos: AAPos
    gene_refstart: NAPos
    gene_refend: NAPos
    gene: GeneText
    curidx: int = -1

    posnas_by_napos: List[
        Tuple[NAPos, List[PosNA]]
    ] = group_posnas_by_napos(posnas)
    size: int = len(posnas_by_napos)
    basepairs: List[BasePair] = []

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

            basepairs.append((
                gene,
                aapos,
                na_and_ins
            ))
    return basepairs


@cython.cfunc
@cython.inline
@cython.returns(list)
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


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def get_comparable_codon(
    codon_posnas: List[List[PosNA]]
) -> Tuple[Tuple[MultiNAText, ...], bool]:
    posnas: List[PosNA]
    codon_nas: Tuple[MultiNAText, ...] = tuple([
        bytes([na for _, _, na, _ in posnas])
        for posnas in codon_posnas
    ])
    is_partial: bool = len(codon_nas) < 3
    return codon_nas, is_partial


@cython.cfunc
@cython.inline
@cython.returns(list)
def group_codons(
    basepairs: List[BasePair]
) -> List[Tuple[GeneText, AAPos, List[List[PosNA]]]]:
    """Group base-pairs into complete codons

    A codon is represented by a nested list. The inner List[PosNA]
    is an individual base-pair with its insertions; the outer
    List[List[PosNA]] is a complete codon
    """
    gene: GeneText
    aapos: AAPos
    prev_gene: GeneText = '\a0FakeGene\a0'
    prev_aapos: AAPos = -1

    bp: List[PosNA]
    codons: List[Tuple[GeneText, AAPos, List[List[PosNA]]]] = []

    for gene, aapos, na_and_ins in basepairs:
        if gene == prev_gene and aapos == prev_aapos:
            codons[-1][2].append(na_and_ins)
        else:
            prev_gene = gene
            prev_aapos = aapos
            codons.append((gene, aapos, [na_and_ins]))
    return codons


@cython.cfunc
@cython.returns(list)
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
    codon: Tuple[MultiNAText, ...]
    is_partial: bool
    totalq: int
    sizeq: int

    genes: List[GeneInterval] = find_overlapped_genes(
        gene_intervals, read_refstart, read_refend)
    basepairs: List[BasePair] = group_basepairs(posnas, genes)

    poscodons: List[PosCodon] = []
    for gene, aapos, codon_posnas in group_codons(basepairs):
        codon, is_partial = get_comparable_codon(codon_posnas)
        if is_partial:
            continue

        totalq = 0
        sizeq = 0
        for pnas in codon_posnas:
            for pna in pnas:
                totalq += pna[3]
                sizeq += 1
        meanq_int = round(totalq / sizeq if totalq else 0)
        if meanq_int < site_quality_cutoff:
            continue

        poscodons.append((gene, aapos, codon, meanq_int))
    return poscodons


@cython.ccall
@cython.returns(list)
def get_poscodons_between(
    samfile: str,
    samfile_start: int,
    samfile_end: int,
    gene_intervals: List[GeneInterval],
    site_quality_cutoff: int = 0
) -> List[Tuple[Header, List[PosCodon]]]:

    read: AlignedSegment
    posnas: List[PosNA]
    poscodons: List[PosCodon]

    results: List[Tuple[Header, List[PosCodon]]] = []

    with pysam.AlignmentFile(samfile, 'rb') as samfp:
        samfp.seek(samfile_start)

        for read in samfp:
            if samfp.tell() >= samfile_end:
                break

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

            results.append((read.query_name, poscodons))
    return results


def iter_poscodons(
    samfile: str,
    genes: List[DerivedFragmentConfig],
    workers: int,
    description: str,
    site_quality_cutoff: int = 0,
    log_format: str = 'text',
    chunk_size: int = 20000,
    **extras: Any
) -> Generator[
    Tuple[Header, List[PosCodon]],
    None,
    None
]:

    total: int
    read: AlignedSegment
    posnas: List[PosNA]
    chunks: List[Tuple[int, int]]
    samfile_begin: int
    samfile_end: int
    one: Tuple[Header, List[PosCodon]]
    poscodons: List[Tuple[Header, List[PosCodon]]]
    pbar: Optional[Union[JsonProgress, tqdm]]

    with pysam.AlignmentFile(samfile, 'rb') as samfp:
        total = samfp.mapped
        if log_format == 'json':
            pbar = JsonProgress(
                total=total, description=description, **extras)
        elif log_format == 'text':
            pbar = tqdm(total=total)
            pbar.set_description('Processing {}'.format(description))

    chunks = chunked_samfile(samfile, chunk_size)
    gene_intervals: List[GeneInterval] = build_gene_intervals(genes)

    with ProcessPoolExecutor(workers) as executor:

        for poscodons in executor.map(
            get_poscodons_between,
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
            if pbar:
                for one in poscodons:
                    yield one
                    pbar.update(1)
            else:
                yield from poscodons
    if pbar:
        pbar.close()
