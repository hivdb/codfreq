from pathlib import Path
from unittest.mock import MagicMock, patch

from typer.testing import CliRunner
import sys

from .mock_postalign import mock_postalign

_POSTALIGN = mock_postalign()
_POSTALIGN.__enter__()
sys.modules.pop("codfreq.codonalign_consensus", None)
sys.modules.pop("codfreq.sam2codfreq", None)
from codfreq.align import app as align_app  # noqa: E402


def test_align_cli_invalid_program(tmp_path: Path) -> None:
    profile = tmp_path / "profile.json"
    profile.write_text("{}", encoding="utf-8")
    runner = CliRunner()
    result = runner.invoke(
        align_app,
        [str(tmp_path), "-p", "invalid", "-r", str(profile)],
    )
    assert result.exit_code != 0


def test_align_cli_invalid_log_format(tmp_path: Path) -> None:
    profile = tmp_path / "profile.json"
    profile.write_text("{}", encoding="utf-8")
    runner = CliRunner()
    result = runner.invoke(
        align_app,
        [
            str(tmp_path),
            "-p",
            "bowtie2",
            "-r",
            str(profile),
            "--log-format",
            "bad",
        ],
    )
    assert result.exit_code != 0


def test_align_cli_enable_profiling(tmp_path: Path) -> None:
    """Profiling flag runs under cProfile and prints stats."""

    profile = tmp_path / "profile.json"
    profile.write_text("{}", encoding="utf-8")
    runner = CliRunner()
    prof_ctx = MagicMock()
    prof_ctx.__enter__.return_value = prof_ctx
    prof_ctx.__exit__.return_value = None
    with (
        patch("cProfile.Profile", return_value=prof_ctx) as prof_cls,
        patch("pstats.Stats") as stats_cls,
        patch("codfreq.align.align") as mock_align,
    ):
        result = runner.invoke(
            align_app,
            [
                str(tmp_path),
                "-p",
                "bowtie2",
                "-r",
                str(profile),
                "--enable-profiling",
            ],
        )
    assert result.exit_code == 0
    prof_cls.assert_called_once()
    mock_align.assert_called_once()
    stats_cls.return_value.print_stats.assert_called_once()
