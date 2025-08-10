from pathlib import Path

from typer.testing import CliRunner

from codfreq.align import app as align_app


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
