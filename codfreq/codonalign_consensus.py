# cython: profile=True

from postalign.utils.group_by_codons import group_by_codons
from postalign.processors.codon_alignment import codon_align
from postalign.models.sequence import Sequence, PositionalSeqStr
from itertools import groupby
from operator import itemgetter


ENCODING = 'UTF-8'
GAP = ord(b'-')
DEL_CODON = bytearray(b'---')
CODON_ALIGN_WINDOW_SIZE = 5


def get_sequence_obj(seq_bytearray, seqid):
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


def assemble_alignment(geneposcdf_lookup, refseq, genedef):
    gene = genedef['geneName']
    gene_refstart = genedef['refStart']
    gene_refend = genedef['refEnd']
    gene_refseq = bytearray()
    gene_queryseq = bytearray()
    refsize = gene_refend - gene_refstart + 1
    first_aa = refsize // 3
    last_aa = 0

    for aapos in range(1, refsize // 3 + 1):
        napos = gene_refstart + aapos * 3 - 4
        ref_codon = refseq[napos:napos + 3]
        geneposcdf = geneposcdf_lookup.get((gene, aapos))
        if geneposcdf:
            cons_codon = bytearray(geneposcdf[0]['codon'], ENCODING)
            first_aa = min(first_aa, aapos)
            last_aa = max(last_aa, aapos)
        else:
            cons_codon = DEL_CODON
        cons_codon_size = len(cons_codon)
        if cons_codon_size < 3:
            cons_codon += bytearray([GAP]) * (3 - cons_codon_size)
        elif cons_codon_size > 3:
            ref_codon += bytearray([GAP]) * (cons_codon_size - 3)
        gene_refseq += ref_codon
        gene_queryseq += cons_codon

    if last_aa == 0:
        return None, None, None, None

    return (get_sequence_obj(gene_refseq, 1),
            get_sequence_obj(gene_queryseq, 2),
            first_aa,
            last_aa)


def codonalign_consensus(codonfreq, ref, genes):
    refseq = bytearray(ref['refSequence'], ENCODING)
    geneposcdf_lookup = {
        genepos: sorted(genecdf, key=lambda cdf: -cdf['count'])
        for genepos, genecdf in
        groupby(codonfreq, key=itemgetter('gene', 'position'))
    }
    for genedef in genes:
        gene = genedef['geneName']

        (gene_refseq_obj,
         gene_queryseq_obj,
         first_aa,
         last_aa) = assemble_alignment(
             geneposcdf_lookup, refseq, genedef
        )

        if last_aa is None:
            continue

        (gene_refseq_obj,
         gene_queryseq_obj) = codon_align(
             gene_refseq_obj, gene_queryseq_obj,
             window_size=CODON_ALIGN_WINDOW_SIZE,
             refstart=first_aa * 3 - 2,
             refend=last_aa * 3
        )
        # if gene == 'S':
        #     print(str(gene_refseq_obj.seqtext))
        #     print(str(gene_queryseq_obj.seqtext))
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
            # if gene == 'S' and aapos in (68, 69, 70):
            #     print(gene, aapos, querycodon)
            geneposcdf[0]['codon'] = str(
                PositionalSeqStr.join(querycodon)
            )
            yield from geneposcdf
