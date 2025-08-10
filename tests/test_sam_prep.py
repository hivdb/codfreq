"""Tests for :mod:`codfreq.sam_prep`."""

from collections import Counter
from typing import Any, Iterator, List

from unittest.mock import MagicMock, patch

from codfreq.sam_prep import count_indel_positions, prepare_sam, squash_gaps


def test_squash_gaps_merges_indels() -> None:
    cig = ((0, 5), (1, 2), (0, 3), (2, 1), (0, 4))
    assert squash_gaps(cig) == ((0, 5), (1, 1), (0, 7))


def test_count_indel_positions() -> None:
    counter: Counter[int] = Counter()
    cig = ((0, 5), (1, 2), (0, 3), (2, 1), (0, 4))
    count_indel_positions(cig, 100, counter)
    assert counter == Counter({105: 1, 108: 1})


def test_squash_gaps_switches_operation() -> None:
    """Insertions longer than nearby deletions flip the operation."""

    cig = ((0, 5), (1, 2), (0, 3), (2, 4), (0, 1))
    assert squash_gaps(cig) == ((0, 5), (2, 2), (0, 4))


def test_squash_gaps_merges_same_operation() -> None:
    """Adjacent insertions merge into a longer insertion."""

    cig = ((0, 5), (1, 2), (0, 3), (1, 1), (0, 4))
    assert squash_gaps(cig) == ((0, 5), (1, 3), (0, 7))


def test_prepare_sam_processes_reads() -> None:
    """Mapped reads receive squashed cigars and are written out."""

    class Read:
        def __init__(
            self,
            unmapped: bool,
            cigar: tuple[tuple[int, int], ...],
            ref_start: int,
        ) -> None:
            self.is_unmapped = unmapped
            self.cigartuples = list(cigar)
            self.reference_start = ref_start

    unmapped = Read(True, (), 0)
    mapped = Read(False, ((0, 5), (1, 2), (0, 3), (2, 1), (0, 4)), 100)

    reads: List[Read] = [unmapped, mapped]
    written: List[Any] = []

    def alignmentfile_factory(*args: Any, **kwargs: Any) -> MagicMock:
        mode = args[1] if len(args) > 1 else kwargs.get("mode", "r")
        af = MagicMock()
        af.__enter__.return_value = af
        af.__exit__.return_value = None
        if mode.startswith("r"):
            af._pos = 0

            def iterate() -> Iterator[Read]:
                for r in reads[af._pos:]:
                    af._pos += 1
                    yield r

            af.__iter__.side_effect = iterate
            af.tell.side_effect = lambda: af._pos
            af.seek.side_effect = lambda pos: setattr(af, "_pos", pos)
        else:
            af.write.side_effect = lambda r: written.append(r)
        return af

    with patch(
        "codfreq.sam_prep.AlignmentFile", side_effect=alignmentfile_factory
    ):
        prepare_sam("in.sam", "out.sam")

    assert len(written) == 2
    assert written[1].cigartuples == [(0, 5), (1, 1), (0, 7)]
