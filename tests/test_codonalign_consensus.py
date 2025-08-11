"""Tests for :mod:`codfreq.codonalign_consensus`."""

from collections import Counter
from unittest.mock import patch
from codfreq.codonalign_consensus import (
    aapos_to_napos,
    assemble_alignment,
    codonalign_consensus,
)


def test_aapos_to_napos_across_ranges() -> None:
    """AA positions map correctly across multiple reference ranges."""
    refranges = [(1, 3), (10, 18)]
    assert aapos_to_napos(1, refranges) == 1
    assert aapos_to_napos(2, refranges) == 10
    assert aapos_to_napos(4, refranges) == 16
    assert aapos_to_napos(5, refranges) == -1


def test_assemble_alignment_handles_codon_sizes() -> None:
    """Consensus assembly pads or trims codons to fit the reference."""
    codonstat = {
        ("F", 1): Counter({b"AA": 1}),
        ("F", 3): Counter({b"AAAA": 1}),
    }
    fragment = {"fragmentName": "F", "refRanges": [(1, 9)]}
    refseq = bytearray(b"AAACCCGGG")

    ref_obj, query_obj, first, last = assemble_alignment(
        codonstat, refseq, fragment
    )
    assert ref_obj is not None and query_obj is not None
    assert ref_obj.seqtext == bytearray(b"AAACCCGGG-")
    assert query_obj.seqtext == bytearray(b"AA----AAAA")
    assert (first, last) == (1, 3)


def test_codonalign_consensus_updates_counters() -> None:
    """Codon alignment config is honored and counters are updated."""
    codonstat = {
        ("frag", 1): Counter({b"AAA": 2}),
        ("skip", 1): Counter({b"CCC": 1}),
    }
    quals = {
        ("frag", 1): Counter({b"AAA": 5}),
        ("skip", 1): Counter({b"CCC": 3}),
    }
    ref = {"fragmentName": "ref", "refSequence": "AAACCC"}
    fragments = [
        {
            "fragmentName": "skip",
            "refRanges": [(1, 3)],
            "codonAlignment": False,
        },
        {
            "fragmentName": "frag",
            "refRanges": [(1, 3)],
            "codonAlignment": [
                {
                    "relRefStart": 0,
                    "relRefEnd": 15,
                    "minGapDistance": 5,
                    "windowSize": 7,
                    "relGapPlacementScore": "0:0-1=1",
                }
            ],
        },
    ]

    codonalign_consensus(codonstat, quals, ref, fragments)

    assert codonstat[("skip", 1)][b"CCC"] == 1
    assert b"AAA" not in codonstat[("frag", 1)]
    assert codonstat[("frag", 1)][b"TTT"] == 2
    assert b"AAA" not in quals[("frag", 1)]
    assert quals[("frag", 1)][b"TTT"] == 5


def test_assemble_alignment_returns_none_without_codons() -> None:
    """Fragments without codons yield ``None`` results."""
    codonstat: dict[tuple[str, int], Counter[bytes]] = {}
    fragment = {"fragmentName": "F", "refRanges": [(1, 3)]}
    refseq = bytearray(b"AAA")

    ref_obj, query_obj, first, last = assemble_alignment(
        codonstat, refseq, fragment
    )

    assert (ref_obj, query_obj, first, last) == (None, None, None, None)


def test_codonalign_consensus_skips_empty_fragment() -> None:
    """Fragments without statistics are ignored."""
    codonstat: dict[tuple[str, int], Counter[bytes]] = {}
    quals: dict[tuple[str, int], Counter[bytes]] = {}
    ref = {"fragmentName": "ref", "refSequence": "AAA"}
    fragments = [{"fragmentName": "frag", "refRanges": [(1, 3)]}]

    result = codonalign_consensus(codonstat, quals, ref, fragments)
    assert result == (codonstat, quals)


def test_codonalign_consensus_alignment_failure_skips_fragment() -> None:
    """Alignment failures lead to early fragment skipping."""
    codonstat = {("frag", 1): Counter({b"AAA": 1})}
    quals = {("frag", 1): Counter({b"AAA": 1})}
    ref = {"fragmentName": "ref", "refSequence": "AAA"}
    fragments = [{"fragmentName": "frag", "refRanges": [(1, 3)]}]

    with patch(
        "codfreq.codonalign_consensus.codon_align",
        return_value=(None, None),
    ):
        codonalign_consensus(codonstat, quals, ref, fragments)


def test_codonalign_consensus_skips_positions_missing_quality() -> None:
    """Codon positions lacking quality scores are ignored."""
    codonstat = {
        ("frag", 1): Counter({b"AAA": 1}),
        ("frag", 2): Counter({b"CCC": 1}),
    }
    quals = {("frag", 1): Counter({b"AAA": 1})}
    ref = {"fragmentName": "ref", "refSequence": "AAA CCC".replace(" ", "")}
    fragments = [{"fragmentName": "frag", "refRanges": [(1, 6)]}]

    codonalign_consensus(codonstat, quals, ref, fragments)
