"""Tests for alignment helper functions."""

from __future__ import annotations

import json
import os
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open

import pytest
import typer

from codfreq.align import (
    find_paired_marker,
    find_paired_fastq_patterns,
    complete_paired_fastqs,
    fastp_preprocess,
    find_paired_fastqs,
    ivar_trim,
    cutadapt_trim,
    align_with_profile,
    align,
    REQUIRED_PROFILE_VERSION,
)
from codfreq.enums import LogFormat, Program  # noqa: E402
from codfreq.codfreq_types import PairedFASTQ, Profile  # noqa: E402
from codfreq.cmdwrappers.fastp import FASTPConfig  # noqa: E402
from codfreq.cmdwrappers.ivar import TrimConfig  # noqa: E402
from codfreq.cmdwrappers.cutadapt import CutadaptConfig  # noqa: E402


def test_find_paired_marker_valid_and_invalid() -> None:
    """Markers are detected only for proper R1/R2 differences."""
    assert find_paired_marker("read_R1.fastq", "read_R2.fastq") == 6
    assert find_paired_marker("read1.fastq", "read3.fastq") == -1


def test_find_paired_marker_rejects_invalid_patterns() -> None:
    """Extra digits or multiple mismatches invalidate the marker."""
    assert find_paired_marker("sample_R11.fastq", "sample_R12.fastq") == -1
    assert find_paired_marker("1x1a", "2x2b") == -1


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


def test_find_paired_fastq_patterns_invalid_pairs() -> None:
    """Pairs with multiple differences are treated as singles."""
    files = ["a_R1.fastq", "b_R2.fastq"]
    patterns = list(find_paired_fastq_patterns(files, autopairing=True))
    assert len(patterns) == 2
    assert {p["pair"][0] for p in patterns} == set(files)
    assert all(p["pair"][1] is None and p["n"] == 1 for p in patterns)


def test_find_paired_fastq_patterns_sorts_pairs() -> None:
    """Input order does not affect output pair ordering."""

    files = ["sample_R2.fastq", "sample_R1.fastq"]
    patterns = list(find_paired_fastq_patterns(files, autopairing=True))
    assert patterns == [
        {
            "name": "sample",
            "pair": ("sample_R1.fastq", "sample_R2.fastq"),
            "n": 2,
        }
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


def test_find_paired_fastqs_scans_directory(tmp_path: Path) -> None:
    """Directories without pairinfo are scanned for FASTQ files."""
    (tmp_path / "sample_R1.fastq").write_text("")
    (tmp_path / "sample_R2.fastq").write_text("")
    (tmp_path / "single.fastq").write_text("")
    results = list(find_paired_fastqs(str(tmp_path), autopairing=True))
    assert len(results) == 2
    pairinfo = tmp_path / "pairinfo.json"
    assert pairinfo.is_file()


def test_ivar_trim_logs_and_calls_wrapper(tmp_path: Path) -> None:
    """ivar.trim is invoked and progress logged in text mode."""
    in_bam = str(tmp_path / "input.bam")
    out_bam = str(tmp_path / "output.bam")
    config: TrimConfig = {}
    with (
        patch("codfreq.align.ivar.trim") as mock_trim,
        patch("codfreq.align.rich.print") as mock_print,
    ):
        ivar_trim(in_bam, out_bam, config, LogFormat.text)
    mock_trim.assert_called_once_with(in_bam, out_bam, **config)
    mock_print.assert_any_call(
        f"Trimming {os.path.basename(in_bam)} using ivar..."
    )


def test_cutadapt_trim_json_logs(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """cutadapt.cutadapt is called and emits JSON progress."""
    merged: PairedFASTQ = {
        "name": "sample",
        "pair": (str(tmp_path / "sample.merged.fastq.gz"), None),
        "n": 1,
    }
    config: CutadaptConfig = {}
    with patch("codfreq.align.cutadapt.cutadapt") as mock_cut:
        result = cutadapt_trim(merged, config, LogFormat.json)
    out = capsys.readouterr().out
    assert '"command": "cutadapt"' in out
    mock_cut.assert_called_once_with(
        merged["pair"][0], result["pair"][0], **config
    )
    assert result["pair"][0].endswith("merged-trimed.fastq.gz")


def test_find_paired_fastq_patterns_chunk_mismatch() -> None:
    """Filename pairs with different chunk counts are treated as singles."""

    files = ["a_extra_R1.fastq", "b_R2.fastq"]
    patterns = list(find_paired_fastq_patterns(files, autopairing=True))
    assert len(patterns) == 2
    assert all(p["pair"][1] is None for p in patterns)


def test_find_paired_fastq_patterns_multiple_marker_diffs() -> None:
    """Pairs with more than one marker difference are rejected."""

    files = ["s_R1_1_x.fastq", "s_R2_2_x.fastq"]
    patterns = list(find_paired_fastq_patterns(files, autopairing=True))
    assert len(patterns) == 2
    assert all(p["pair"][1] is None for p in patterns)


def test_ivar_trim_json_logs(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """ivar.trim emits structured JSON progress when requested."""

    in_bam = str(tmp_path / "in.bam")
    out_bam = str(tmp_path / "out.bam")
    config: TrimConfig = {}
    with patch("codfreq.align.ivar.trim") as mock_trim:
        ivar_trim(in_bam, out_bam, config, LogFormat.json)
    out = capsys.readouterr().out
    assert '"command": "ivar"' in out
    mock_trim.assert_called_once_with(in_bam, out_bam, **config)


def test_cutadapt_trim_text_logs(tmp_path: Path) -> None:
    """cutadapt.cutadapt emits human-readable progress in text mode."""

    merged: PairedFASTQ = {
        "name": "sample",
        "pair": (str(tmp_path / "sample.merged.fastq.gz"), None),
        "n": 1,
    }
    config: CutadaptConfig = {}
    with (
        patch("codfreq.align.cutadapt.cutadapt") as mock_cut,
        patch("codfreq.align.rich.print") as mock_print,
    ):
        result = cutadapt_trim(merged, config, LogFormat.text)
    mock_print.assert_any_call("Trimming sample using cutadapt...")
    mock_cut.assert_called_once_with(
        merged["pair"][0], result["pair"][0], **config
    )
    assert result["pair"][0].endswith("merged-trimed.fastq.gz")


def test_align_with_profile_replaces_without_trim(tmp_path: Path) -> None:
    """When no trimming is configured files are renamed after alignment."""

    paired: PairedFASTQ = {"name": "samp", "pair": ("r1.fq", "r2.fq"), "n": 2}
    profile = Profile.model_validate(
        {
            "version": "1",
            "fragmentConfig": [{"fragmentName": "F", "refSequence": "AAA"}],
            "sequenceAssemblyConfig": [],
        }
    )
    with (
        patch("codfreq.align.fastp_preprocess", return_value=paired),
        patch("codfreq.align.get_refinit", return_value=lambda x: None),
        patch(
            "codfreq.align.get_align",
            return_value=MagicMock(),
        ) as get_align,
        patch("codfreq.align.os.replace") as mock_replace,
        patch("codfreq.align.rich.print"),
    ):
        align_with_profile(
            paired,
            Program.minimap2,
            profile,
            LogFormat.text,
            fastp_config={},
            cutadapt_config=None,
            ivar_trim_config=None,
        )
    get_align.return_value.assert_called_once()
    assert mock_replace.call_count == 3


def test_align_with_profile_trims_and_logs_json(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Alignment path uses cutadapt and ivar trimming when configured."""

    paired: PairedFASTQ = {"name": "samp", "pair": ("r1.fq", "r2.fq"), "n": 2}
    profile = Profile.model_validate(
        {
            "version": "1",
            "fragmentConfig": [{"fragmentName": "F", "refSequence": "AAA"}],
            "sequenceAssemblyConfig": [],
        }
    )
    with (
        patch("codfreq.align.fastp_preprocess", return_value=paired),
        patch("codfreq.align.cutadapt_trim", return_value=paired) as mock_cut,
        patch("codfreq.align.get_refinit", return_value=lambda x: None),
        patch("codfreq.align.get_align", return_value=MagicMock()),
        patch("codfreq.align.ivar_trim") as mock_ivar,
    ):
        align_with_profile(
            paired,
            Program.bowtie2,
            profile,
            LogFormat.json,
            fastp_config={},
            cutadapt_config={},
            ivar_trim_config={},
        )
    out = capsys.readouterr().out
    assert '"op": "alignment"' in out
    mock_cut.assert_called_once()
    mock_ivar.assert_called_once()


def test_align_aborts_on_profile_version_mismatch(tmp_path: Path) -> None:
    """Profiles with wrong version trigger a ``typer.Abort``."""

    profile = tmp_path / "p.json"
    profile.write_text("{}")
    with (
        profile.open() as fp,
        patch("codfreq.align.json.load", return_value={"version": "old"}),
        patch("codfreq.align.rich.print") as mock_print,
    ):
        with pytest.raises(typer.Abort):
            align(tmp_path, Program.bowtie2, fp, 1, LogFormat.text, True)
    mock_print.assert_called_once()


def test_align_runs_pipeline(tmp_path: Path) -> None:
    """Core pipeline orchestrates alignment and codfreq generation."""

    profile = tmp_path / "p.json"
    profile.write_text("{}")
    profile_obj = {
        "version": REQUIRED_PROFILE_VERSION,
        "fragmentConfig": [],
        "sequenceAssemblyConfig": [],
    }
    pairobj = {"name": "samp", "pair": ("r1", "r2"), "n": 2}
    with (
        profile.open() as fp,
        patch("codfreq.align.json.load", return_value=profile_obj),
        patch("codfreq.align.find_paired_fastqs", return_value=[pairobj]),
        patch("codfreq.align.fastp.load_config", return_value={}),
        patch("codfreq.align.cutadapt.load_config", return_value=None),
        patch("codfreq.align.ivar.load_trim_config", return_value=None),
        patch("codfreq.align.align_with_profile") as mock_align_profile,
        patch(
            "codfreq.align.name_codfreq",
            return_value=str(tmp_path / "out.codfreq"),
        ),
        patch("codfreq.align.open", mock_open(), create=True),
        patch("codfreq.align.csv.DictWriter") as mock_writer,
        patch(
            "codfreq.align.sam2codfreq_all",
            return_value=[{"codon": b"AAA"}],
        ),
        patch(
            "codfreq.align.create_untrans_region_consensus"
        ) as mock_consensus,
    ):
        align(tmp_path, Program.minimap2, fp, 1, LogFormat.text, True)
    mock_align_profile.assert_called_once()
    mock_writer.return_value.writeheader.assert_called_once()
    mock_writer.return_value.writerow.assert_called_once()
    mock_consensus.assert_called_once()
