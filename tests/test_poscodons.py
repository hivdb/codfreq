"""Tests for :mod:`codfreq.poscodons`."""

from codfreq.poscodons import (
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
