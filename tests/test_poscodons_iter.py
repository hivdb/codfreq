"""Tests for :func:`codfreq.poscodons.iter_poscodons`."""

import sys
import types
import importlib
from array import array
from typing import Any, Iterator

# Stub pysam module
pysam_stub = types.ModuleType("pysam")


class AlignedSegment:
    """Minimal aligned read stub."""

    def __init__(self, name: str, seq: str, qual: list[int]) -> None:
        self.query_name = name
        self.query_sequence = seq
        self.query_qualities = array("B", qual)
        self.reference_start = 0
        self.reference_end = len(seq)

    def get_aligned_pairs(self, _with_seq: bool) -> list[tuple[int, int]]:
        return [(i, i) for i in range(len(self.query_sequence))]


class AlignmentFile:
    """Iterator over preset reads supporting ``seek`` and ``tell``."""

    reads: list[AlignedSegment] = []

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self._pos = 0

    def __enter__(self) -> "AlignmentFile":
        return self

    def __exit__(self, *exc: Any) -> None:  # pragma: no cover - no-op
        return None

    def __iter__(self) -> Iterator[AlignedSegment]:
        while self._pos < len(self.reads):
            yield self.reads[self._pos]
            self._pos += 1

    def seek(self, pos: int) -> None:
        self._pos = pos

    def tell(self) -> int:
        return self._pos


pysam_stub.AlignmentFile = AlignmentFile  # type: ignore[attr-defined]
pysam_stub.AlignedSegment = AlignedSegment  # type: ignore[attr-defined]
sys.modules["pysam"] = pysam_stub

# Stub cython decorators
cython_stub = types.ModuleType("cython")


def _decorator(*dargs: Any, **dkwargs: Any) -> Any:  # pragma: no cover
    if dargs and callable(dargs[0]):
        return dargs[0]
    return lambda func: func


cython_stub.cfunc = _decorator  # type: ignore[attr-defined]
cython_stub.inline = _decorator  # type: ignore[attr-defined]
cython_stub.returns = lambda *a, **k: _decorator  # type: ignore[attr-defined]
sys.modules["cython"] = cython_stub

import codfreq.poscodons as poscodons  # noqa: E402
importlib.reload(poscodons)
from codfreq.codfreq_types import FragmentInterval  # noqa: E402
from codfreq.poscodons import iter_poscodons  # noqa: E402


def test_iter_poscodons_yields_codons() -> None:
    """Reads are converted to positional codons."""

    AlignmentFile.reads = [AlignedSegment("r1", "ATGATT", [30] * 6)]
    frags: list[FragmentInterval] = [([(1, 6)], "frag")]
    result = list(iter_poscodons("sample.sam", 0, 1, frags))
    assert result == [
        (
            "r1",
            [
                ("frag", 1, b"ATG", 30),
                ("frag", 2, b"ATT", 30),
            ],
        )
    ]


def test_iter_poscodons_stops_when_past_end() -> None:
    """Iteration halts once the end position is exceeded."""

    AlignmentFile.reads = [
        AlignedSegment("r1", "AAA", [30, 30, 30]),
        AlignedSegment("r2", "GGG", [30, 30, 30]),
    ]
    frags: list[FragmentInterval] = [([(1, 3)], "frag")]
    result = list(iter_poscodons("sample.sam", 0, 0, frags))
    assert result == [("r1", [("frag", 1, b"AAA", 30)])]


def test_iter_poscodons_skips_empty_reads() -> None:
    """Reads lacking sequence are ignored."""

    AlignmentFile.reads = [
        AlignedSegment("r1", "", []),
        AlignedSegment("r2", "TTT", [30, 30, 30]),
    ]
    frags: list[FragmentInterval] = [([(1, 3)], "frag")]
    result = list(iter_poscodons("sample.sam", 0, 5, frags))
    assert result == [("r2", [("frag", 1, b"TTT", 30)])]
