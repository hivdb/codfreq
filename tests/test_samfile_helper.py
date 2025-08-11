"""Tests for :func:`codfreq.samfile_helper.chunked_samfile`."""

import importlib
import sys
import types
from typing import Any, Iterator
from unittest.mock import patch
import pytest

pysam_stub = types.ModuleType("pysam")


class AlignmentFile:  # pragma: no cover - stub
    """Minimal ``pysam.AlignmentFile`` replacement."""

    n_records = 0

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self.pos = 0
        self.n_records = AlignmentFile.n_records

    def __enter__(self) -> "AlignmentFile":
        return self

    def __exit__(self, *exc: Any) -> None:
        return None

    def __iter__(self) -> Iterator[object]:
        for _ in range(self.n_records):
            self.pos += 1
            yield object()

    def tell(self) -> int:
        return self.pos


pysam_stub.AlignmentFile = AlignmentFile  # type: ignore[attr-defined]

cython_stub = types.ModuleType("cython")


def _decorator(*args: Any, **kwargs: Any) -> Any:  # pragma: no cover - no-op
    if args and callable(args[0]):
        return args[0]
    return lambda func: func


cython_stub.ccall = _decorator  # type: ignore[attr-defined]
cython_stub.inline = _decorator  # type: ignore[attr-defined]
cython_stub.returns = lambda *a, **k: _decorator  # type: ignore[attr-defined]

_patch = patch.dict(sys.modules, {"pysam": pysam_stub, "cython": cython_stub})
_patch.start()

import codfreq.samfile_helper as samfile_helper  # noqa: E402
importlib.reload(samfile_helper)
from codfreq.samfile_helper import chunked_samfile  # noqa: E402


@pytest.fixture(scope="module", autouse=True)
def _cleanup_modules() -> Iterator[None]:
    """Remove stubbed modules after tests.

    Yields:
        Iterator[None]: ``None``.
    """
    yield
    _patch.stop()


def test_chunked_samfile_groups_offsets() -> None:
    """Offsets are grouped according to the chunk size."""

    AlignmentFile.n_records = 12
    chunks = chunked_samfile("sample.sam", chunk_size=5)
    assert chunks == [(0, 5), (5, 10), (10, 12)]
