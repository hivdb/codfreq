import sys
import types
import importlib
from typing import Any, Iterator

pysam_stub = types.ModuleType("pysam")


class AlignmentFile:  # pragma: no cover - stub
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
sys.modules["pysam"] = pysam_stub

cython_stub = types.ModuleType("cython")


def _decorator(*args: Any, **kwargs: Any) -> Any:  # pragma: no cover - no-op
    if args and callable(args[0]):
        return args[0]
    return lambda func: func


cython_stub.ccall = _decorator  # type: ignore[attr-defined]
cython_stub.inline = _decorator  # type: ignore[attr-defined]
cython_stub.returns = lambda *a, **k: _decorator  # type: ignore[attr-defined]
sys.modules["cython"] = cython_stub

import codfreq.samfile_helper as samfile_helper  # noqa: E402
importlib.reload(samfile_helper)
from codfreq.samfile_helper import chunked_samfile  # noqa: E402


def test_chunked_samfile_groups_offsets() -> None:
    AlignmentFile.n_records = 12
    chunks = chunked_samfile("sample.sam", chunk_size=5)
    assert chunks == [(0, 5), (5, 10), (10, 12)]
