"""Tests for :mod:`codfreq.codonutils`."""

from codfreq.codonutils import expand_ambiguous_na, translate_codon


def test_expand_ambiguous_na() -> None:
    """Ambiguous bases expand to all possible nucleotides."""

    assert expand_ambiguous_na(ord(b"R")) == b"AG"
    assert expand_ambiguous_na(ord(b"A")) == b"A"


def test_translate_codon_with_ambiguity() -> None:
    """``translate_codon`` resolves ambiguous NA symbols."""

    assert translate_codon(b"ATG") == b"M"
    assert translate_codon(b"ATN") == b"IM"
    # cached lookup returns the same result
    assert translate_codon(b"ATN") == b"IM"
