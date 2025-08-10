import cython  # type: ignore
import pysam  # type: ignore
from pysam import AlignedSegment  # type: ignore
from collections.abc import Generator
from collections import defaultdict

from .codfreq_types import (
    FragmentInterval,
    Header,
    NAChar,
    AAPos,
    NAPos,
    NAPosRange,
    CodonText
)
from .posnas import iter_single_read_posnas, PosNA

#                                          Qual
#                                           v
PosCodon = tuple[Header, AAPos, CodonText, int]
BasePair = tuple[AAPos, list[PosNA]]


@cython.cfunc
@cython.inline
@cython.returns(list)
def group_posnas_by_napos(
    posnas: list[PosNA]
) -> list[tuple[NAPos, list[PosNA]]]:
    """Group PosNA entries by nucleotide position."""

    prev_pos: int = -1
    by_napos: list[tuple[NAPos, list[PosNA]]] = []
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
    posnas: list[PosNA],
    fragment_intervals: list[FragmentInterval]
) -> list[tuple[Header, list[BasePair]]]:
    """Group positional nucleotides into codon base pairs."""

    napos: NAPos
    na_and_ins: list[PosNA]
    aapos: AAPos
    frag_refranges: list[NAPosRange]
    fragment_name: Header

    posnas_by_napos: list[
        tuple[NAPos, list[PosNA]]
    ] = group_posnas_by_napos(posnas)
    basepairs: defaultdict[Header, list[BasePair]] = defaultdict(list)

    for napos, na_and_ins in posnas_by_napos:
        for frag_refranges, fragment_name in fragment_intervals:
            rel_napos0 = 0
            for start, end in frag_refranges:
                if napos >= start and napos <= end:
                    aapos = (rel_napos0 + napos - start) // 3 + 1
                    basepairs[fragment_name].append((
                        aapos,
                        na_and_ins
                    ))
                rel_napos0 += end - start + 1
    return list(basepairs.items())


@cython.cfunc
@cython.inline
@cython.returns(list)
def find_intersected_fragments(
    fragment_intervals: list[FragmentInterval],
    read_refstart: NAPos,
    read_refend: NAPos
) -> list[FragmentInterval]:
    frag_refranges: list[NAPosRange]
    fragment_name: Header
    filtered: list[FragmentInterval] = []
    for frag_refranges, fragment_name in fragment_intervals:
        if all(read_refend < start for start, _ in frag_refranges):
            continue
        if all(read_refstart > end for _, end in frag_refranges):
            continue
        filtered.append((frag_refranges, fragment_name))
    return filtered


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def get_comparable_codon(
    codon_posnas: list[list[PosNA]]
) -> tuple[CodonText, bool]:
    posnas: list[PosNA]
    codon_chars: list[NAChar] = []
    num_bps: int = 0

    for posnas in codon_posnas:
        num_bps += 1
        for posna in posnas:
            codon_chars.append(posna[2])

    is_partial: bool = num_bps < 3
    return bytes(codon_chars), is_partial


@cython.cfunc
@cython.inline
@cython.returns(list)
def group_codons(
    basepairs: list[tuple[Header, list[BasePair]]]
) -> list[tuple[Header, AAPos, list[list[PosNA]]]]:
    """Group base-pairs into complete codons

    A codon is represented by a nested list. The inner List[PosNA]
    is an individual base-pair with its insertions; the outer
    List[List[PosNA]] is a complete codon
    """
    aapos: AAPos
    fragment_name: Header
    fragment_bps: list[BasePair]
    codons: list[tuple[Header, AAPos, list[list[PosNA]]]] = []

    for fragment_name, fragment_bps in basepairs:
        prev_aapos: AAPos = -1
        for aapos, na_and_ins in fragment_bps:
            if aapos == prev_aapos:
                codons[-1][2].append(na_and_ins)
            else:
                prev_aapos = aapos
                codons.append((fragment_name, aapos, [na_and_ins]))
    return codons


@cython.cfunc
@cython.inline
@cython.returns(list)
def posnas2poscodons(
    posnas: list[PosNA],
    fragment_intervals: list[FragmentInterval],
    read_refstart: int,  # 1-based first aligned refpos
    read_refend: int,    # 1-based last aligned refpos
    site_quality_cutoff: int
) -> list[PosCodon]:
    """Convert positional nucleotides to codons with quality scores."""

    meanq_int: int
    fragment_name: Header
    aapos: AAPos
    codon_posnas: list[list[PosNA]]
    codon: CodonText
    is_partial: bool
    totalq: int
    sizeq: int

    fragments: list[FragmentInterval] = find_intersected_fragments(
        fragment_intervals, read_refstart, read_refend)
    basepairs: list[
        tuple[Header, list[BasePair]]
    ] = group_basepairs(posnas, fragments)

    poscodons: list[PosCodon] = []
    for fragment_name, aapos, codon_posnas in group_codons(basepairs):
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

        poscodons.append((fragment_name, aapos, codon, meanq_int))
    return poscodons


def iter_poscodons(
    samfile: str,
    samfile_start: int,
    samfile_end: int,
    fragment_intervals: list[FragmentInterval],
    site_quality_cutoff: int = 0
) -> Generator[tuple[Header | None, list[PosCodon]], None, None]:
    """Retrieve poscodons from given SAM/BAM file position range"""

    read: AlignedSegment
    posnas: list[PosNA]
    poscodons: list[PosCodon]

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

            poscodons = posnas2poscodons(
                posnas,
                fragment_intervals,
                read.reference_start + 1,  # pysam has 0-based numbering
                read.reference_end,  # "reference_end points to one past the
                                     #  last aligned residue."
                site_quality_cutoff
            )

            yield read.query_name, poscodons
