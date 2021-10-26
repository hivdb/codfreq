from typing import Dict, List, Set
from .codfreq_types import (
    CodonText,
    AAChar, MultiAAText,
    NAChar, MultiNAText
)


CODON_TABLE: Dict[CodonText, MultiAAText] = {
    b'TTT': b'F',
    b'TTC': b'F',
    b'TTA': b'L',
    b'TTG': b'L',

    b'CTT': b'L',
    b'CTC': b'L',
    b'CTA': b'L',
    b'CTG': b'L',

    b'ATT': b'I',
    b'ATC': b'I',
    b'ATA': b'I',
    b'ATG': b'M',

    b'GTT': b'V',
    b'GTC': b'V',
    b'GTA': b'V',
    b'GTG': b'V',

    b'TCT': b'S',
    b'TCC': b'S',
    b'TCA': b'S',
    b'TCG': b'S',

    b'CCT': b'P',
    b'CCC': b'P',
    b'CCA': b'P',
    b'CCG': b'P',

    b'ACT': b'T',
    b'ACC': b'T',
    b'ACA': b'T',
    b'ACG': b'T',

    b'GCT': b'A',
    b'GCC': b'A',
    b'GCA': b'A',
    b'GCG': b'A',

    b'TAT': b'Y',
    b'TAC': b'Y',

    b'CAT': b'H',
    b'CAC': b'H',
    b'CAA': b'Q',
    b'CAG': b'Q',

    b'AAT': b'N',
    b'AAC': b'N',
    b'AAA': b'K',
    b'AAG': b'K',

    b'GAT': b'D',
    b'GAC': b'D',
    b'GAA': b'E',
    b'GAG': b'E',

    b'TGT': b'C',
    b'TGC': b'C',
    b'TGG': b'W',

    b'CGT': b'R',
    b'CGC': b'R',
    b'CGA': b'R',
    b'CGG': b'R',

    b'AGT': b'S',
    b'AGC': b'S',
    b'AGA': b'R',
    b'AGG': b'R',

    b'GGT': b'G',
    b'GGC': b'G',
    b'GGA': b'G',
    b'GGG': b'G',

    b'TAA': b'*',
    b'TGA': b'*',
    b'TAG': b'*',
}

REVERSE_CODON_TABLE: Dict[AAChar, List[CodonText]] = {}
for codon, aa in CODON_TABLE.items():
    REVERSE_CODON_TABLE.setdefault(aa[0], []).append(codon)


AMBIGUOUS_NAS: Dict[NAChar, MultiNAText] = {
    ord(b'W'): b'AT',
    ord(b'S'): b'CG',
    ord(b'M'): b'AC',
    ord(b'K'): b'GT',
    ord(b'R'): b'AG',
    ord(b'Y'): b'CT',
    ord(b'B'): b'CGT',
    ord(b'D'): b'AGT',
    ord(b'H'): b'ACT',
    ord(b'V'): b'ACG',
    ord(b'N'): b'ACGT'
}


def expand_ambiguous_na(na: NAChar) -> MultiNAText:
    return AMBIGUOUS_NAS.get(na, bytes([na]))


def translate_codon(nas: MultiNAText) -> MultiAAText:
    aas_text: bytes
    nas = nas.replace(b'-', b'N')[:3]
    if nas in CODON_TABLE:
        return CODON_TABLE[nas]
    aas: Set[int] = set()
    for na0 in AMBIGUOUS_NAS.get(nas[0], nas[0:1]):
        for na1 in AMBIGUOUS_NAS.get(nas[1], nas[1:2]):
            for na2 in AMBIGUOUS_NAS.get(nas[2], nas[2:3]):
                aas |= set(
                    CODON_TABLE[bytes([na0, na1, na2])]
                )
    CODON_TABLE[nas] = aas_text = bytes((sorted(aas)))
    return aas_text
