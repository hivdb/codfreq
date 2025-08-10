"""Tests for alignment helper functions."""

from __future__ import annotations

import json
import os
from pathlib import Path
from unittest.mock import patch
import sys

import pytest

from .mock_postalign import mock_postalign

_POSTALIGN = mock_postalign()
_POSTALIGN.__enter__()
sys.modules.pop("codfreq.codonalign_consensus", None)
sys.modules.pop("codfreq.sam2codfreq", None)
from codfreq.align import (  # noqa: E402
    find_paired_marker,
    find_paired_fastq_patterns,
    complete_paired_fastqs,
    fastp_preprocess,
    find_paired_fastqs,
)
from codfreq.enums import LogFormat  # noqa: E402
from codfreq.codfreq_types import PairedFASTQ  # noqa: E402
from codfreq.cmdwrappers.fastp import FASTPConfig  # noqa: E402


def test_find_paired_marker_valid_and_invalid() -> None:
    """Markers are detected only for proper R1/R2 differences."""
    assert find_paired_marker("read_R1.fastq", "read_R2.fastq") == 6
    assert find_paired_marker("read1.fastq", "read3.fastq") == -1


def test_find_paired_fastq_patterns_autopairing() -> None:
    """Autopairing groups matching pairs and singles."""
    files = ["sample_R1.fastq", "sample_R2.fastq", "single.fastq"]
    patterns = list(find_paired_fastq_patterns(files, autopairing=True))
    assert patterns == [
        {
            "name": "sample",
            "pair": ("sample_R1.fastq", "sample_R2.fastq"),
            "n": 2,
        },
        {"name": "single", "pair": ("single.fastq", None), "n": 1},
    ]


def test_find_paired_fastq_patterns_no_autopair() -> None:
    """Without autopairing every file is treated as single-ended."""
    files = ["a.fastq", "b.fastq"]
    patterns = list(find_paired_fastq_patterns(files, autopairing=False))
    assert patterns == [
        {"name": "a", "pair": ("a.fastq", None), "n": 1},
        {"name": "b", "pair": ("b.fastq", None), "n": 1},
    ]


def test_complete_paired_fastqs_expands_paths() -> None:
    """Relative paths are joined with the directory path."""
    pairs: list[PairedFASTQ] = [
        {"name": "a", "pair": ("r1.fq", "r2.fq"), "n": 2}
    ]
    result = list(complete_paired_fastqs(pairs, "/work"))
    assert result == [
        {
            "name": "/work/a",
            "pair": ("/work/r1.fq", "/work/r2.fq"),
            "n": 2,
        }
    ]


def test_fastp_preprocess_invokes_wrapper(tmp_path: Path) -> None:
    """Fastp wrapper is called and merged filename returned."""
    pfq: PairedFASTQ = {
        "name": "sample",
        "pair": (str(tmp_path / "r1.fq"), str(tmp_path / "r2.fq")),
        "n": 2,
    }
    config: FASTPConfig = {}
    with (
        patch("codfreq.align.fastp.fastp") as mock_fastp,
        patch("codfreq.align.rich.print") as mock_print,
    ):
        result = fastp_preprocess(pfq, config, LogFormat.text)
    mock_fastp.assert_called_once()
    mock_print.assert_any_call("Pre-processing sample using fastp...")
    assert result["pair"][0].endswith("sample.merged.fastq.gz")


def test_fastp_preprocess_json_logs(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """JSON log format prints structured progress messages."""
    pfq: PairedFASTQ = {
        "name": "sample",
        "pair": (str(tmp_path / "r1.fq"), str(tmp_path / "r2.fq")),
        "n": 2,
    }
    config: FASTPConfig = {}
    with patch("codfreq.align.fastp.fastp") as mock_fastp:
        result = fastp_preprocess(pfq, config, LogFormat.json)
    out = capsys.readouterr().out
    assert '"op": "preprocess"' in out
    assert result["pair"][0].endswith("sample.merged.fastq.gz")
    mock_fastp.assert_called_once()


def test_find_paired_fastqs_reads_pairinfo(tmp_path: Path) -> None:
    """Existing pairinfo file is loaded and expanded."""
    data = [{"name": "samp", "pair": ["a.fastq", None], "n": 1}]
    pairinfo = tmp_path / "pairinfo.json"
    pairinfo.write_text(json.dumps(data))
    results = list(find_paired_fastqs(str(tmp_path), autopairing=False))
    assert results == [
        {
            "name": os.path.join(str(tmp_path), "samp"),
            "pair": (os.path.join(str(tmp_path), "a.fastq"), None),
            "n": 1,
        }
    ]
