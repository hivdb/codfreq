from more_itertools import chunked  # type: ignore

from .codfreq_types import RefAAs, Sequence, SeqText, NAPos
from .fastareader import load
from .codonutils import translate_codon


def get_refaas(reference: SeqText, reference_start: NAPos = 1) -> RefAAs:
    refseq: Sequence
    refaas: RefAAs = {}
    with open(reference) as fp:
        refseq, = load(fp)
        refseq_text: bytes = refseq['sequence']
        refseq_text = refseq_text[reference_start - 1:]
        for pos0, codon in enumerate(chunked(refseq_text, 3)):
            refaas[pos0 + 1] = translate_codon(bytes(codon))
    return refaas
