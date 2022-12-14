import cython  # type: ignore
from postalign.utils.group_by_codons import group_by_codons
from postalign.processors.codon_alignment import (
    codon_align,
    parse_gap_placement_score
)
from postalign.models.sequence import (
    Sequence,
    NAPosition
)
from operator import itemgetter
from typing import (
    Optional,
    List,
    Dict,
    Tuple,
    Iterator,
    Union,
    Literal,
    Counter
)

from .sam2codfreq_types import CodonCounterByFragPos
from .codfreq_types import (
    AAPos,
    NAPos,
    NAPosRange,
    CodonText,
    Header,
    CodFreqRow,
    MainFragmentConfig,
    DerivedFragmentConfig,
    CodonAlignmentConfig
)


ENCODING = 'UTF-8'
GAP = ord(b'-')
DEL_CODON = b'---'
CODON_ALIGN_WINDOW_SIZE = 10
CODON_ALIGN_MIN_GAP_DISTANCE = 30


def get_sequence_obj(
    seq_bytearray: bytearray,
    seqid: int
) -> Sequence:
    return Sequence(
        header='seq{}'.format(seqid),
        description='',
        seqtext=NAPosition.init_from_bytes(seq_bytearray),
        seqid=seqid,
        seqtype=NAPosition,
        abs_seqstart=0,
        skip_invalid=True
    )


@cython.cfunc
@cython.inline
@cython.returns(int)
def aapos_to_napos(
    aapos: AAPos,
    refranges: List[NAPosRange]
) -> NAPos:
    max_rel_napos = 0
    for start, end in refranges:
        max_rel_napos += end - start + 1
        max_aapos = max_rel_napos // 3
        if aapos <= max_aapos:
            rel_napos = aapos * 3 - 2
            napos_offset = max_rel_napos - rel_napos
            return end - napos_offset
    return -1


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def assemble_alignment(
    codonstat_by_fragpos: CodonCounterByFragPos,
    refseq: bytearray,
    fragment: DerivedFragmentConfig
) -> Tuple[
    Optional[Sequence],
    Optional[Sequence],
    Optional[AAPos],
    Optional[AAPos]
]:
    aapos: AAPos
    napos: NAPos
    ref_codon: bytearray
    cons_codon: bytearray
    cons_codon_bytes: bytes
    cons_codon_size: int
    codons: Optional[Counter[CodonText]]
    fragment_name: Header = fragment['fragmentName']
    frag_refranges: List[NAPosRange] = fragment['refRanges']
    frag_refseq: bytearray = bytearray()
    frag_queryseq: bytearray = bytearray()
    refsize: int = sum(end - start + 1 for start, end in frag_refranges)
    first_aa: AAPos = refsize // 3
    last_aa: AAPos = 0

    for aapos in range(1, refsize // 3 + 1):
        napos = aapos_to_napos(aapos, frag_refranges)
        if napos == -1:
            # aapos is out of range, should we raise error?
            continue
        ref_codon = refseq[napos - 1:napos + 2]
        codons = codonstat_by_fragpos.get((fragment_name, aapos))
        if codons:
            ((cons_codon_bytes, _),) = codons.most_common(1)
            first_aa = min(first_aa, aapos)
            last_aa = max(last_aa, aapos)
        else:
            cons_codon_bytes = DEL_CODON
        cons_codon = bytearray(cons_codon_bytes)
        cons_codon_size = len(cons_codon_bytes)
        if cons_codon_size < 3:
            cons_codon.extend([GAP] * (3 - cons_codon_size))
        elif cons_codon_size > 3:
            ref_codon.extend([GAP] * (cons_codon_size - 3))
        frag_refseq.extend(ref_codon)
        frag_queryseq.extend(cons_codon)

    if last_aa == 0:
        return None, None, None, None

    return (
        get_sequence_obj(frag_refseq, 1),
        get_sequence_obj(frag_queryseq, 2),
        first_aa,
        last_aa
    )


genepos_getter = itemgetter('gene', 'position')
codon_getter = itemgetter('codon')
order_getter = itemgetter('count', 'total_quality_score')


@cython.ccall
@cython.inline
@cython.returns(tuple)
def codonalign_consensus(
    codonstat_by_fragpos: CodonCounterByFragPos,
    qualities_by_fragpos: CodonCounterByFragPos,
    ref: MainFragmentConfig,
    fragments: List[DerivedFragmentConfig],
) -> Tuple[CodonCounterByFragPos, CodonCounterByFragPos]:
    fragment_name: Header
    fragment: DerivedFragmentConfig
    frag_refseq_obj: Optional[Sequence]
    frag_queryseq_obj: Optional[Sequence]
    first_aa: Optional[AAPos]
    last_aa: Optional[AAPos]
    seq_ref_start: NAPos
    aapos0: AAPos
    aapos: AAPos
    refstart: NAPos
    refend: NAPos
    refcodon: List[NAPosition]
    querycodon: List[NAPosition]
    codons: Optional[Counter[CodonText]]
    cdfs: Iterator[CodFreqRow]
    cdf_list: List[CodFreqRow]
    codon: CodonText

    refseq: bytearray = bytearray(ref['refSequence'], ENCODING)
    for fragment in fragments:
        fragment_name = fragment['fragmentName']
        codon_align_config: Optional[
            Union[Literal[False], List[CodonAlignmentConfig]]
        ] = fragment.get('codonAlignment')

        if codon_align_config is False:
            # skip this gene if explicitly defined codonAlignment=False
            continue

        if not codon_align_config:
            codon_align_config = [{}]

        # assemble consensus codon reads into pairwise alignment
        (frag_refseq_obj,
         frag_queryseq_obj,
         first_aa,
         last_aa) = assemble_alignment(
             codonstat_by_fragpos, refseq, fragment
        )

        if (
            frag_refseq_obj is None or
            frag_queryseq_obj is None or
            first_aa is None or
            last_aa is None
        ):
            continue

        # convert AA positions to NA positions
        seq_refstart = first_aa * 3 - 2
        seq_refend = last_aa * 3

        # apply codon alignment (CDA) to pairwise alignment
        for cda_config in codon_align_config:
            refstart = cda_config.get('relRefStart', seq_refstart)
            refend = cda_config.get('relRefEnd', seq_refend)

            # Codon alignment shouldn't exceed query sequence boundary
            if refstart < seq_refend and refend > seq_refstart:
                refstart = max(refstart, seq_refstart)
                refend = min(refend, seq_refend)

            # Load minGapDistance, windowSize and gapPlacementScore from config
            min_gap_distance = cda_config.get(
                'minGapDistance'
            ) or CODON_ALIGN_MIN_GAP_DISTANCE
            window_size = cda_config.get(
                'windowSize'
            ) or CODON_ALIGN_WINDOW_SIZE
            gap_placement_score: Dict[
                int, Dict[Tuple[int, int], int]
            ] = parse_gap_placement_score(
                cda_config.get('relGapPlacementScore') or '')

            # perform postalign's codon_align
            (frag_refseq_obj,
             frag_queryseq_obj) = codon_align(
                 frag_refseq_obj,
                 frag_queryseq_obj,
                 min_gap_distance=min_gap_distance,
                 window_size=window_size,
                 gap_placement_score=gap_placement_score,
                 ref_start=refstart,
                 ref_end=refend
            )

        if (
            frag_refseq_obj is None or
            frag_queryseq_obj is None
        ):
            continue

        for aapos0, (refcodon, querycodon) in enumerate(zip(
            *group_by_codons(
                frag_refseq_obj.seqtext,
                frag_queryseq_obj.seqtext
            )
        )):
            aapos = aapos0 + 1
            oldcodon: CodonText
            newcodon: CodonText = NAPosition.as_bytes(querycodon)

            codons = codonstat_by_fragpos.get((fragment_name, aapos))
            quas = qualities_by_fragpos.get((fragment_name, aapos))
            if not codons or not quas:
                continue

            # replace the consensus codon to codon aligned codon
            ((oldcodon, _),) = codons.most_common(1)
            if newcodon != oldcodon:
                codons[newcodon] += codons.pop(oldcodon)
                quas[newcodon] += quas.pop(oldcodon)

    return codonstat_by_fragpos, qualities_by_fragpos
