"""Tests for :mod:`codfreq.posnas`."""

from array import array

from codfreq.posnas import GAP, iter_single_read_posnas


def test_iter_single_read_posnas_basic() -> None:
    """Simple alignments yield sequential positions."""

    seq = "AC"
    qual = array("B", [10, 20])
    pairs = [(0, 0), (1, 1)]
    assert iter_single_read_posnas(seq, qual, pairs) == [
        (1, 0, ord("A"), 10),
        (2, 0, ord("C"), 20),
    ]


def test_iter_single_read_posnas_insertion_deletion() -> None:
    """Insertions get an index and deletions emit the gap char."""

    seq = "ACG"
    qual = array("B", [10, 20, 30])
    pairs = [(0, 0), (1, None), (2, 1), (None, 2)]
    posnas = iter_single_read_posnas(seq, qual, pairs)
    assert posnas == [
        (1, 0, ord("A"), 10),
        (1, 1, ord("C"), 20),
        (2, 0, ord("G"), 30),
        (3, 0, GAP, 30),
    ]


def test_iter_single_read_posnas_trims_trailing_insertions() -> None:
    """Trailing insertions are removed from the output."""

    seq = "ACG"
    qual = array("B", [10, 20, 30])
    pairs = [(0, 0), (1, 1), (2, None)]
    assert iter_single_read_posnas(seq, qual, pairs) == [
        (1, 0, ord("A"), 10),
        (2, 0, ord("C"), 20),
    ]


def test_iter_single_read_posnas_skips_leading_insertions() -> None:
    """Insertions before the first reference base are ignored."""

    seq = "AC"
    qual = array("B", [10, 20])
    pairs = [(0, None), (1, 0)]
    assert iter_single_read_posnas(seq, qual, pairs) == [
        (1, 0, ord("C"), 20),
    ]


def test_iter_single_read_posnas_trims_multiple_trailing_insertions() -> None:
    """Multiple trailing insertions are removed from the output."""

    seq = "ACGT"
    qual = array("B", [10, 20, 30, 40])
    pairs = [(0, 0), (1, 1), (2, None), (3, None)]
    assert iter_single_read_posnas(seq, qual, pairs) == [
        (1, 0, ord("A"), 10),
        (2, 0, ord("C"), 20),
    ]


def test_iter_single_read_posnas_defaults_quality() -> None:
    """Quality scores default to ``1`` when absent."""

    seq = "AC"
    pairs = [(0, 0), (1, 1)]
    assert iter_single_read_posnas(seq, None, pairs) == [
        (1, 0, ord("A"), 1),
        (2, 0, ord("C"), 1),
    ]
