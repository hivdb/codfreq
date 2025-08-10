"""Tests for codon-alignment helpers and consensus assembly."""

import sys
import types
from collections import Counter
from typing import Any

# Stub out postalign dependencies used by codonalign_consensus
postalign = types.ModuleType("postalign")
postalign.__path__ = []
utils = types.ModuleType("postalign.utils")
utils.__path__ = []
processors = types.ModuleType("postalign.processors")
processors.__path__ = []
models = types.ModuleType("postalign.models")
models.__path__ = []

sys.modules["postalign"] = postalign
sys.modules["postalign.utils"] = utils
sys.modules["postalign.processors"] = processors
sys.modules["postalign.models"] = models

# group_by_codons stub
utils_group = types.ModuleType("postalign.utils.group_by_codons")


def group_by_codons(
    seq1: bytearray,
    seq2: bytearray,
) -> tuple[list[bytearray], list[bytearray]]:
    """Group two sequences into codon-sized chunks."""

    def chunk(seq: bytearray) -> list[bytearray]:
        return [seq[i:i + 3] for i in range(0, len(seq), 3)]

    return chunk(seq1), chunk(seq2)


utils_group.group_by_codons = group_by_codons  # type: ignore[attr-defined]
sys.modules["postalign.utils.group_by_codons"] = utils_group

# codon_alignment stubs
processors_codon = types.ModuleType("postalign.processors.codon_alignment")


def codon_align(refseq: Any, queryseq: Any, **_kwargs: Any) -> tuple[Any, Any]:
    """Pretend to realign codons by mutating the query sequence."""
    if len(queryseq.seqtext) >= 3:
        queryseq.seqtext[:3] = b"TTT"
    return refseq, queryseq


def parse_gap_placement_score(_text: str) -> dict:
    """Return an empty gap-placement map for tests."""
    return {}


processors_codon.codon_align = codon_align  # type: ignore[attr-defined]
processors_codon.parse_gap_placement_score = (  # type: ignore[attr-defined]
    parse_gap_placement_score
)
sys.modules["postalign.processors.codon_alignment"] = processors_codon

# sequence model stubs
models_seq = types.ModuleType("postalign.models.sequence")


class NAPosition:
    """Minimal representation of a nucleotide position."""

    @staticmethod
    def init_from_bytes(b: bytearray) -> bytearray:
        return bytearray(b)

    @staticmethod
    def as_bytes(seq: bytearray) -> bytes:
        return bytes(seq)


class Sequence:
    """Simple sequence container mimicking postalign's model."""

    def __init__(
        self,
        header: str,
        description: str,
        seqtext: bytearray,
        seqid: int,
        seqtype: Any,
        abs_seqstart: int,
        skip_invalid: bool,
    ) -> None:
        self.header = header
        self.description = description
        self.seqtext = seqtext
        self.seqid = seqid
        self.seqtype = seqtype
        self.abs_seqstart = abs_seqstart
        self.skip_invalid = skip_invalid


models_seq.NAPosition = NAPosition  # type: ignore[attr-defined]
models_seq.Sequence = Sequence  # type: ignore[attr-defined]
sys.modules["postalign.models.sequence"] = models_seq

sys.modules.pop("codfreq.codonalign_consensus", None)

from codfreq.codonalign_consensus import (  # noqa: E402
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
