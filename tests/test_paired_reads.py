"""Tests for :mod:`codfreq.paired_reads`."""

from pathlib import Path

from unittest.mock import MagicMock, patch

from codfreq.paired_reads import iter_paired_reads


def test_iter_paired_reads_groups_by_query_names(tmp_path: Path) -> None:
    """Reads sharing a query name are grouped together."""

    reads = [
        MagicMock(query_name="r1"),
        MagicMock(query_name="r1"),
        MagicMock(query_name=None),
        MagicMock(query_name="r2"),
    ]

    alignment_mock = MagicMock()
    alignment_mock.__enter__.return_value = alignment_mock
    alignment_mock.__exit__.return_value = None
    alignment_mock.fetch.return_value = reads

    with patch(
        "codfreq.paired_reads.pysam.AlignmentFile",
        return_value=alignment_mock,
    ):
        samfile = tmp_path / "reads.sam"
        samfile.write_text("", encoding="utf-8")
        pairs = iter_paired_reads(str(samfile))

    assert [name for name, _ in pairs] == ["r1", "r2"]
    assert len(pairs[0][1]) == 2
