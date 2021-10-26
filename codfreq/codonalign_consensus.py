from postalign.utils.group_by_codons import group_by_codons  # type: ignore
from postalign.processors.codon_alignment import codon_align  # type: ignore
from postalign.models.sequence import (  # type: ignore
    Sequence,
    PositionalSeqStr
)
from itertools import groupby
from operator import itemgetter
from typing import Generator, Optional, List, Dict, Tuple, Iterator

from .codfreq_types import (
    AAPos,
    NAPos,
    CodonText,
    GeneText,
    CodFreqRow,
    MainFragmentConfig,
    DerivedFragmentConfig
)


ENCODING = 'UTF-8'
GAP = ord(b'-')
DEL_CODON = b'---'
CODON_ALIGN_WINDOW_SIZE = 5


def get_sequence_obj(
    seq_bytearray: bytearray,
    seqid: int
) -> Sequence:
    return Sequence(
        header='seq{}'.format(seqid),
        description='',
        seqtext=PositionalSeqStr.init_from_nastring(
            bytes(seq_bytearray).decode(ENCODING)
        ),
        seqid=seqid,
        seqtype='NA',
        abs_seqstart=0,
        skip_invalid=True
    )


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

    return (get_sequence_obj(gene_refseq, 1),
            get_sequence_obj(gene_queryseq, 2),
            first_aa,
            last_aa)


genepos_getter = itemgetter('gene', 'position')
codon_getter = itemgetter('codon')
count_getter = itemgetter('count')


def codonalign_consensus(
    codonfreq: Generator[CodFreqRow, None, None],
    ref: MainFragmentConfig,
    genes: List[DerivedFragmentConfig],
) -> Generator[CodFreqRow, None, None]:
    gene: GeneText
    genedef: DerivedFragmentConfig
    gene_ref_seq_obj: Optional[Sequence]
    gene_queryseq_obj: Optional[Sequence]
    first_aa: Optional[AAPos]
    last_aa: Optional[AAPos]
    aapos0: AAPos
    aapos: AAPos
    refcodon: List[PositionalSeqStr]
    querycodon: List[PositionalSeqStr]
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
        genepos: sorted(genecdf, key=count_getter, reverse=True)
        for genepos, genecdf in
        groupby(codonfreq, key=genepos_getter)
    }
    for genedef in genes:
        gene = genedef['geneName']

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

        (gene_refseq_obj,
         gene_queryseq_obj) = codon_align(
             gene_refseq_obj, gene_queryseq_obj,
             window_size=CODON_ALIGN_WINDOW_SIZE,
             refstart=first_aa * 3 - 2,
             refend=last_aa * 3
        )
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
            geneposcdf[0]['codon'] = str(
                PositionalSeqStr.join(querycodon)
            ).encode(ENCODING)

            # merge same codons
            merged_geneposcdf = []
            for codon, cdfs in groupby(
                sorted(geneposcdf, key=codon_getter),
                codon_getter
            ):
                cdf_list = list(cdfs)
                cdf_list[0]['count'] = sum(cdf['count'] for cdf in cdf_list)
                cdf_list[0]['total_quality_score'] = sum(
                    cdf['total_quality_score'] for cdf in cdf_list
                )
                merged_geneposcdf.append(cdf_list[0])
            yield from sorted(
                merged_geneposcdf,
                key=count_getter,
                reverse=True
            )
