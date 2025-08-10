"""Tests for I/O helpers in :mod:`codfreq.posnas`."""

import sys
import types
import importlib
from array import array
from typing import Any, Iterator
from unittest.mock import patch, MagicMock

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
        self.mapped = len(self.reads)

    def __enter__(self) -> "AlignmentFile":
        return self

    def __exit__(self, *exc: Any) -> None:  # pragma: no cover - no-op
        return None

    def __iter__(self) -> Iterator[AlignedSegment]:
        while self._pos < len(self.reads):
            read = self.reads[self._pos]
            self._pos += 1
            yield read

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
    iter_posnas,
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


def test_get_posnas_between_skips_empty_reads() -> None:
    """Reads without sequence are ignored."""

    AlignmentFile.reads = [
        AlignedSegment("r1", "AC", [10, 20], [(0, 0), (1, 1)]),
        AlignedSegment("r2", "", [], []),
    ]
    result = get_posnas_between("sample.sam", 0, 2)
    assert result == [("r1", [(1, 0, ord("A"), 10), (2, 0, ord("C"), 20)])]


def test_get_posnas_in_genome_region_skips_empty_reads() -> None:
    """Reads without sequence are ignored."""

    AlignmentFile.reads = [
        AlignedSegment("r1", "AC", [10, 20], [(0, 0), (1, 1)]),
        AlignedSegment("r2", "", [], []),
    ]
    result = get_posnas_in_genome_region("sample.sam", "ref", 1, 5)
    assert result == [("r1", [(1, 0, ord("A"), 10), (2, 0, ord("C"), 20)])]


def test_get_posnas_in_genome_region_filters_quality() -> None:
    """Low-quality bases are removed from genomic region results."""

    AlignmentFile.reads = [
        AlignedSegment("r1", "AC", [10, 20], [(0, 0), (1, 1)]),
    ]
    result = get_posnas_in_genome_region(
        "sample.sam", "ref", 1, 5, site_quality_cutoff=15
    )
    assert result == [("r1", [(2, 0, ord("C"), 20)])]


def test_iter_posnas_reports_progress() -> None:
    """Chunks are processed and progress bar updated per read."""

    AlignmentFile.reads = [
        AlignedSegment("r1", "AC", [10, 20], [(0, 0), (1, 1)]),
        AlignedSegment("r2", "GT", [30, 40], [(0, 2), (1, 3)]),
    ]

    class DummyExecutor:
        def __init__(self, _workers: int) -> None:
            pass

        def __enter__(self) -> "DummyExecutor":
            return self

        def __exit__(self, *exc: Any) -> None:
            return None

        def map(self, func: Any, *iterables: Any) -> Any:
            return map(func, *iterables)

    with (
        patch("codfreq.posnas.ProcessPoolExecutor", DummyExecutor),
        patch("codfreq.posnas.JsonProgress") as MockProgress,
        patch("codfreq.posnas.chunked_samfile", return_value=[(0, 1), (1, 2)]),
    ):
        results = list(
            iter_posnas(
                "sample.sam",
                workers=1,
                description="reads",
                log_format="json",
                chunk_size=1,
            )
        )
    mock_bar = MockProgress.return_value
    assert results == [
        ("r1", [(1, 0, ord("A"), 10), (2, 0, ord("C"), 20)]),
        ("r2", [(3, 0, ord("G"), 30), (4, 0, ord("T"), 40)]),
    ]
    assert mock_bar.update.call_count == 2
    mock_bar.close.assert_called_once()


def test_iter_posnas_without_progress() -> None:
    """Results are yielded directly when no progress bar is configured."""

    AlignmentFile.reads = [
        AlignedSegment("r1", "AC", [10, 20], [(0, 0), (1, 1)]),
    ]

    class DummyExecutor:
        def __init__(self, _workers: int) -> None:
            pass

        def __enter__(self) -> "DummyExecutor":
            return self

        def __exit__(self, *exc: Any) -> None:
            return None

        def map(self, func: Any, *iterables: Any) -> Any:
            return map(func, *iterables)

    with (
        patch("codfreq.posnas.ProcessPoolExecutor", DummyExecutor),
        patch("codfreq.posnas.chunked_samfile", return_value=[(0, 1)]),
    ):
        results = list(
            iter_posnas(
                "sample.sam",
                workers=1,
                description="reads",
                log_format="quiet",
                chunk_size=1,
            )
        )
    assert results == [
        ("r1", [(1, 0, ord("A"), 10), (2, 0, ord("C"), 20)])
    ]


def test_iter_posnas_text_progress() -> None:
    """Text progress uses ``tqdm`` and sets the description."""
    AlignmentFile.reads = [
        AlignedSegment("r1", "AC", [10, 20], [(0, 0), (1, 1)]),
        AlignedSegment("r2", "GT", [30, 40], [(0, 2), (1, 3)]),
    ]

    class DummyExecutor:
        def __init__(self, _workers: int) -> None:
            pass

        def __enter__(self) -> "DummyExecutor":
            return self

        def __exit__(self, *exc: Any) -> None:
            return None

        def map(self, func: Any, *iterables: Any) -> Any:
            return map(func, *iterables)

    bar = MagicMock()
    with (
        patch("codfreq.posnas.ProcessPoolExecutor", DummyExecutor),
        patch("codfreq.posnas.tqdm", return_value=bar) as mock_tqdm,
        patch("codfreq.posnas.chunked_samfile", return_value=[(0, 1), (1, 2)]),
    ):
        results = list(
            iter_posnas(
                "sample.sam",
                workers=1,
                description="reads",
                log_format="text",
                chunk_size=1,
            )
        )
    mock_tqdm.assert_called_once_with(total=2)
    bar.set_description.assert_called_once_with("Processing reads")
    assert bar.update.call_count == 2
    bar.close.assert_called_once()
    assert results == [
        ("r1", [(1, 0, ord("A"), 10), (2, 0, ord("C"), 20)]),
        ("r2", [(3, 0, ord("G"), 30), (4, 0, ord("T"), 40)]),
    ]
