import cython  # type: ignore
import pysam  # type: ignore
from pysam import AlignedSegment  # type: ignore
from typing import List, Dict, Tuple, Generator, Optional
from collections import defaultdict

from .codfreq_types import (
    FragmentInterval,
    Header,
    NAChar,
    AAPos,
    NAPos,
    CodonText
)
from .posnas import iter_single_read_posnas, PosNA

#                                          Qual
#                                           v
PosCodon = Tuple[Header, AAPos, CodonText, int]
BasePair = Tuple[AAPos, List[PosNA]]


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
    fragment_intervals: List[FragmentInterval]
) -> List[Tuple[Header, List[BasePair]]]:
    """Group same base-pair posnas (NA and ins) by its fragment AA position"""

    napos: NAPos
    na_and_ins: List[PosNA]
    aapos: AAPos
    frag_refstart: NAPos
    frag_refend: NAPos
    fragment_name: Header

    posnas_by_napos: List[
        Tuple[NAPos, List[PosNA]]
    ] = group_posnas_by_napos(posnas)
    basepairs: Dict[Header, List[BasePair]] = defaultdict(list)

    for napos, na_and_ins in posnas_by_napos:
        for frag_refstart, frag_refend, fragment_name in fragment_intervals:
            if napos < frag_refstart or napos > frag_refend:
                continue
            aapos = (napos - frag_refstart) // 3 + 1
            basepairs[fragment_name].append((
                aapos,
                na_and_ins
            ))
    return list(basepairs.items())


@cython.cfunc
@cython.inline
@cython.returns(list)
def find_intersected_fragments(
    fragment_intervals: List[FragmentInterval],
    read_refstart: NAPos,
    read_refend: NAPos
) -> List[FragmentInterval]:
    frag_refstart: NAPos
    frag_refend: NAPos
    fragment_name: Header
    filtered: List[FragmentInterval] = []
    for frag_refstart, frag_refend, fragment_name in fragment_intervals:
        if read_refend < frag_refstart:
            break
        if read_refstart > frag_refend:
            continue
        filtered.append((frag_refstart, frag_refend, fragment_name))
    return filtered


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def get_comparable_codon(
    codon_posnas: List[List[PosNA]]
) -> Tuple[CodonText, bool]:
    posnas: List[PosNA]
    codon_chars: List[NAChar] = []
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
    basepairs: List[Tuple[Header, List[BasePair]]]
) -> List[Tuple[Header, AAPos, List[List[PosNA]]]]:
    """Group base-pairs into complete codons

    A codon is represented by a nested list. The inner List[PosNA]
    is an individual base-pair with its insertions; the outer
    List[List[PosNA]] is a complete codon
    """
    aapos: AAPos
    fragment_name: Header
    fragment_bps: List[BasePair]

    bp: List[PosNA]
    codons: List[Tuple[Header, AAPos, List[List[PosNA]]]] = []

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
    posnas: List[PosNA],
    fragment_intervals: List[FragmentInterval],
    read_refstart: int,  # 1-based first aligned refpos
    read_refend: int,    # 1-based last aligned refpos
    site_quality_cutoff: int
) -> List[PosCodon]:
    meanq: List[int]
    meanq_int: int
    fragment_name: Header
    aapos: AAPos
    codon_posnas: List[List[PosNA]]
    codon: CodonText
    is_partial: bool
    totalq: int
    sizeq: int

    fragments: List[FragmentInterval] = find_intersected_fragments(
        fragment_intervals, read_refstart, read_refend)
    basepairs: List[
        Tuple[Header, List[BasePair]]
    ] = group_basepairs(posnas, fragments)

    poscodons: List[PosCodon] = []
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
    fragment_intervals: List[FragmentInterval],
    site_quality_cutoff: int = 0
) -> Generator[Tuple[Optional[Header], List[PosCodon]], None, None]:
    """Retrieve poscodons from given SAM/BAM file position range"""

    read: AlignedSegment
    posnas: List[PosNA]
    poscodons: List[PosCodon]

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
