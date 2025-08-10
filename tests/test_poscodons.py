"""Tests for :mod:`codfreq.poscodons`."""

from codfreq.poscodons import (
    find_intersected_fragments,
    group_posnas_by_napos,
    posnas2poscodons,
)
from codfreq.codfreq_types import FragmentInterval


def test_group_posnas_by_napos() -> None:
    """PosNAs are grouped by their reference positions."""

    posnas = [(1, 0, ord("A"), 10), (1, 1, ord("C"), 20), (2, 0, ord("G"), 30)]
    grouped = group_posnas_by_napos(posnas)
    assert grouped == [
        (1, [(1, 0, ord("A"), 10), (1, 1, ord("C"), 20)]),
        (2, [(2, 0, ord("G"), 30)]),
    ]


def test_posnas2poscodons_basic() -> None:
    """Positional nucleotides convert to codons with mean quality."""

    posnas = [
        (1, 0, ord("A"), 30),
        (2, 0, ord("T"), 30),
        (3, 0, ord("G"), 30),
        (4, 0, ord("A"), 20),
        (5, 0, ord("T"), 20),
        (6, 0, ord("T"), 20),
    ]
    frags: list[FragmentInterval] = [([(1, 6)], "frag")]
    poscodons = posnas2poscodons(posnas, frags, 1, 6, 0)
    assert poscodons == [
        ("frag", 1, b"ATG", 30),
        ("frag", 2, b"ATT", 20),
    ]


def test_find_intersected_fragments_skips_outside() -> None:
    """Fragments entirely before or after a read are ignored."""

    frags: list[FragmentInterval] = [([(10, 20)], "frag")]
    assert find_intersected_fragments(frags, 1, 5) == []
    assert find_intersected_fragments(frags, 21, 30) == []


def test_posnas2poscodons_filters_partial_and_low_quality() -> None:
    """Incomplete or low-quality codons are excluded."""

    frags: list[FragmentInterval] = [([(1, 3)], "frag")]
    partial = [(1, 0, ord("A"), 30), (2, 0, ord("T"), 30)]
    assert posnas2poscodons(partial, frags, 1, 2, 0) == []

    low_quality = [
        (1, 0, ord("A"), 10),
        (2, 0, ord("T"), 10),
        (3, 0, ord("G"), 10),
    ]
    assert posnas2poscodons(low_quality, frags, 1, 3, 15) == []
