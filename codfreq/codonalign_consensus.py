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
from itertools import groupby
from operator import itemgetter
from typing import Optional, List, Dict, Tuple, Iterator

from .codfreq_types import (
    AAPos,
    NAPos,
    CodonText,
    GeneText,
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
@cython.returns(tuple)
def assemble_alignment(
    geneposcdf_lookup: Dict[
        Tuple[GeneText, AAPos],
        List[CodFreqRow]
    ],
    refseq: bytearray,
    genedef: DerivedFragmentConfig
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
    cons_codon_size: int
    geneposcdf: Optional[List[CodFreqRow]]
    gene: GeneText = genedef['geneName']
    gene_refstart: NAPos = genedef['refStart']
    gene_refend: NAPos = genedef['refEnd']
    gene_refseq: bytearray = bytearray()
    gene_queryseq: bytearray = bytearray()
    refsize: int = gene_refend - gene_refstart + 1
    first_aa: AAPos = refsize // 3
    last_aa: AAPos = 0

    for aapos in range(1, refsize // 3 + 1):
        napos = gene_refstart + aapos * 3 - 4
        ref_codon = refseq[napos:napos + 3]
        geneposcdf = geneposcdf_lookup.get((gene, aapos))
        if geneposcdf:
            cons_codon = bytearray(geneposcdf[0]['codon'])
            first_aa = min(first_aa, aapos)
            last_aa = max(last_aa, aapos)
        else:
            cons_codon = bytearray(DEL_CODON)
        cons_codon_size = len(cons_codon)
        if cons_codon_size < 3:
            cons_codon.extend([GAP] * (3 - cons_codon_size))
        elif cons_codon_size > 3:
            ref_codon.extend([GAP] * (cons_codon_size - 3))
        gene_refseq.extend(ref_codon)
        gene_queryseq.extend(cons_codon)

    if last_aa == 0:
        return None, None, None, None

    return (
        get_sequence_obj(gene_refseq, 1),
        get_sequence_obj(gene_queryseq, 2),
        first_aa,
        last_aa
    )


genepos_getter = itemgetter('gene', 'position')
codon_getter = itemgetter('codon')
order_getter = itemgetter('count', 'total_quality_score')


@cython.ccall
@cython.inline
@cython.returns(list)
def codonalign_consensus(
    codonfreq: List[CodFreqRow],
    ref: MainFragmentConfig,
    genes: List[DerivedFragmentConfig],
) -> List[CodFreqRow]:
    gene: GeneText
    genedef: DerivedFragmentConfig
    gene_refseq_obj: Optional[Sequence]
    gene_queryseq_obj: Optional[Sequence]
    first_aa: Optional[AAPos]
    last_aa: Optional[AAPos]
    seq_ref_start: NAPos
    aapos0: AAPos
    aapos: AAPos
    refstart: NAPos
    refend: NAPos
    ref_offset: NAPos
    refcodon: List[NAPosition]
    querycodon: List[NAPosition]
    geneposcdf: Optional[List[CodFreqRow]]
    merged_geneposcdf: List[CodFreqRow]
    cdfs: Iterator[CodFreqRow]
    cdf_list: List[CodFreqRow]
    codon: CodonText
    refseq: bytearray = bytearray(ref['refSequence'], ENCODING)
    geneposcdf_lookup: Dict[
        Tuple[GeneText, AAPos],
        List[CodFreqRow]
    ] = {
        genepos: sorted(genecdf, key=order_getter, reverse=True)
        for genepos, genecdf in
        groupby(codonfreq, key=genepos_getter)
    }
    rows: List[CodFreqRow] = []
    for genedef in genes:
        gene = genedef['geneName']
        codon_align_config: Optional[
            List[CodonAlignmentConfig]
        ] = genedef.get('codonAlignment')
        if not codon_align_config:
            codon_align_config = [{
                'refStart': genedef['refStart'],
                'refEnd': genedef['refEnd']
            }]

        (gene_refseq_obj,
         gene_queryseq_obj,
         first_aa,
         last_aa) = assemble_alignment(
             geneposcdf_lookup, refseq, genedef
        )

        if (
            gene_refseq_obj is None or
            gene_queryseq_obj is None or
            first_aa is None or
            last_aa is None
        ):
            continue

        seq_refstart = first_aa * 3 - 2
        seq_refend = last_aa * 3
        ref_offset = genedef['refStart'] - 1

        for cda_config in codon_align_config:
            refstart = cda_config['refStart'] - ref_offset
            refend = cda_config['refEnd'] - ref_offset
            if refstart < seq_refend and refend > seq_refstart:
                refstart = max(refstart, seq_refstart)
                refend = min(refend, seq_refend)
            min_gap_distance = cda_config.get(
                'minGapDistance'
            ) or CODON_ALIGN_MIN_GAP_DISTANCE
            window_size = cda_config.get(
                'windowSize'
            ) or CODON_ALIGN_WINDOW_SIZE
            gap_placement_score: Dict[
                int, Dict[Tuple[int, int], int]
            ] = parse_gap_placement_score(
                cda_config.get('gapPlacementScore') or '')
            # position in gap_placement_score should all
            # be subtracted by ref_offset
            gap_placement_score = {
                indel: {
                    (pos - ref_offset, size): score
                    for (pos, size), score in scores.items()
                }
                for indel, scores in gap_placement_score.items()
            }
            (gene_refseq_obj,
             gene_queryseq_obj) = codon_align(
                 gene_refseq_obj, gene_queryseq_obj,
                 min_gap_distance=min_gap_distance,
                 window_size=window_size,
                 gap_placement_score=gap_placement_score,
                 refstart=refstart,
                 refend=refend,
                 check_boundary=False
            )

        if (
            gene_refseq_obj is None or
            gene_queryseq_obj is None
        ):
            continue

        for aapos0, (refcodon, querycodon) in enumerate(zip(
            *group_by_codons(
                gene_refseq_obj.seqtext,
                gene_queryseq_obj.seqtext
            )
        )):
            aapos = aapos0 + 1
            geneposcdf = geneposcdf_lookup.get((gene, aapos))
            if not geneposcdf:
                continue
            geneposcdf[0]['codon'] = NAPosition.as_bytes(querycodon)

            # merge same codons
            merged_geneposcdf = []
            for codon, cdfs in groupby(
                sorted(geneposcdf, key=codon_getter),
                codon_getter
            ):
                cdf_list = list(cdfs)
                cdf_list[0]['count'] = sum(
                    [cdf['count'] for cdf in cdf_list]
                )
                cdf_list[0]['total_quality_score'] = sum(
                    [cdf['total_quality_score'] for cdf in cdf_list]
                )
                merged_geneposcdf.append(cdf_list[0])
            rows.extend(sorted(
                merged_geneposcdf,
                key=order_getter,
                reverse=True
            ))
    return rows
