from postalign.processors.codon_alignment import codon_align
# from postalign.processors.apply_frameshift import (
#     apply_frameshift_to_single_seq, parse_frameshift
# )

CODON_ALIGN_WINDOW_SIZE = 5


def codonalign_caller(refseq, queryseq, refconfig):
    config = refconfig.get('codonAlignments', [])
    for refstart, refend in config:
        refseq, queryseq = codon_align(
            refseq, queryseq,
            window_size=CODON_ALIGN_WINDOW_SIZE,
            refstart=refstart,
            refend=refend
        )
    return refseq, queryseq
