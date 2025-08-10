"""Tests for I/O helpers in :mod:`codfreq.posnas`."""

import sys
import types
import importlib
from array import array
from typing import Any, Iterator

# Stub pysam module
pysam_stub = types.ModuleType("pysam")


class AlignedSegment:
    """Minimal aligned read stub."""

    def __init__(
        self,
        name: str,
        seq: str,
        qual: list[int],
        pairs: list[tuple[int | None, int | None]],
    ) -> None:
        self.query_name = name
        self.query_sequence = seq
        self.query_qualities = array("B", qual)
        self._pairs = pairs

    def get_aligned_pairs(
        self, _with_seq: bool
    ) -> list[tuple[int | None, int | None]]:
        return self._pairs


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

    def fetch(self, *args: Any) -> Iterator[AlignedSegment]:
        yield from self.reads


pysam_stub.AlignmentFile = AlignmentFile  # type: ignore[attr-defined]
pysam_stub.AlignedSegment = AlignedSegment  # type: ignore[attr-defined]
sys.modules["pysam"] = pysam_stub

# Stub cython decorators
cython_stub = types.ModuleType("cython")


def _decorator(*dargs: Any, **dkwargs: Any) -> Any:  # pragma: no cover
    if dargs and callable(dargs[0]):
        return dargs[0]
    return lambda func: func


cython_stub.ccall = _decorator  # type: ignore[attr-defined]
cython_stub.inline = _decorator  # type: ignore[attr-defined]
cython_stub.returns = lambda *a, **k: _decorator  # type: ignore[attr-defined]
sys.modules["cython"] = cython_stub

import codfreq.posnas as posnas  # noqa: E402
importlib.reload(posnas)
from codfreq.posnas import (  # noqa: E402
    get_posnas_between,
    get_posnas_in_genome_region,
)


def test_get_posnas_between_filters_quality() -> None:
    """Low-quality bases are filtered out."""

    AlignmentFile.reads = [
        AlignedSegment("r1", "AC", [10, 20], [(0, 0), (1, 1)]),
        AlignedSegment("r2", "GT", [5, 30], [(0, 2), (1, 3)]),
    ]
    result = get_posnas_between("sample.sam", 0, 2, site_quality_cutoff=15)
    assert result == [
        ("r1", [(2, 0, ord("C"), 20)]),
        ("r2", [(4, 0, ord("T"), 30)]),
    ]


def test_get_posnas_in_genome_region_skips_empty_reads() -> None:
    """Reads without sequence are ignored."""

    AlignmentFile.reads = [
        AlignedSegment("r1", "AC", [10, 20], [(0, 0), (1, 1)]),
        AlignedSegment("r2", "", [], []),
    ]
    result = get_posnas_in_genome_region("sample.sam", "ref", 1, 5)
    assert result == [("r1", [(1, 0, ord("A"), 10), (2, 0, ord("C"), 20)])]
