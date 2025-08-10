"""Tests for :mod:`codfreq.paired_reads`."""

import sys
import types
from pathlib import Path

import pytest


def test_iter_paired_reads_groups_by_query_names(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reads sharing a query name are grouped together."""

    pysam_stub = types.ModuleType("pysam")

    class AlignedSegment:
        def __init__(self, name: str | None) -> None:
            self.query_name = name

    class AlignmentFile:
        def __init__(self, *args: object, **kwargs: object) -> None:
            pass

        def __enter__(self) -> "AlignmentFile":  # noqa: D401
            return self

        def __exit__(self, *exc: object) -> None:  # noqa: D401
            return None

        def fetch(self) -> list[AlignedSegment]:  # noqa: D401
            return [
                AlignedSegment("r1"),
                AlignedSegment("r1"),
                AlignedSegment(None),
                AlignedSegment("r2"),
            ]

    pysam_stub.AlignmentFile = AlignmentFile  # type: ignore[attr-defined]
    pysam_stub.AlignedSegment = AlignedSegment  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "pysam", pysam_stub)

    from codfreq.paired_reads import iter_paired_reads

    samfile = tmp_path / "reads.sam"
    samfile.write_text("", encoding="utf-8")
    pairs = iter_paired_reads(str(samfile))
    assert [name for name, _ in pairs] == ["r1", "r2"]
    assert len(pairs[0][1]) == 2
