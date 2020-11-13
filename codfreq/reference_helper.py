from more_itertools import chunked

from . import fastareader
from .codonutils import translate_codon


def get_refaas(reference, reference_start=1):
    refaas = {}
    with open(reference) as fp:
        refseq, = fastareader.load(fp)
        refseq = refseq['sequence']
        refseq = refseq[reference_start - 1:]
        for pos0, codon in enumerate(chunked(refseq, 3)):
            refaas[pos0 + 1] = translate_codon(''.join(codon))
    return refaas
