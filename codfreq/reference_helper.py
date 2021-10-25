from more_itertools import chunked

from .codfreq_types import RefAAs, Sequence
from .fastareader import load
from .codonutils import translate_codon


def get_refaas(reference: str, reference_start: int = 1) -> RefAAs:
    refseq: Sequence
    refaas: RefAAs = {}
    with open(reference) as fp:
        refseq, = load(fp)
        refseq_text: str = refseq['sequence']
        refseq_text = refseq_text[reference_start - 1:]
        for pos0, codon in enumerate(chunked(refseq, 3)):
            refaas[pos0 + 1] = translate_codon(''.join(codon))
    return refaas
