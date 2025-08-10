from pathlib import Path
import gzip

from unittest.mock import patch

from typer.testing import CliRunner

from codfreq.compress_codfreq import app
from codfreq.cmdwrappers import pigz


def test_compress_codfreq_cli(tmp_path: Path) -> None:
    cf = tmp_path / "sample.codfreq"
    cf.write_text("content", encoding="utf-8")
    ut = tmp_path / "sample.untrans.json"
    ut.write_text("[]", encoding="utf-8")
    with patch.object(pigz, "compress", side_effect=gzip.compress):
        runner = CliRunner()
        result = runner.invoke(app, [str(tmp_path)])
    assert result.exit_code == 0
    assert (tmp_path / "sample.codfreq.gz").exists()


def test_compress_cli_json_logging(tmp_path: Path) -> None:
    """The CLI logs JSON when requested."""

    cf = tmp_path / "sample.codfreq"
    cf.write_text("content", encoding="utf-8")
    ut = tmp_path / "sample.untrans.json"
    ut.write_text("[]", encoding="utf-8")
    with patch.object(pigz, "compress", side_effect=gzip.compress):
        runner = CliRunner()
        result = runner.invoke(app, [str(tmp_path), "--log-format", "json"])
    assert result.exit_code == 0
    assert "compress-codfreq" in result.stdout
    assert (tmp_path / "sample.codfreq.gz").exists()


def test_compress_codfreq_cli_invalid_log_format(tmp_path: Path) -> None:
    runner = CliRunner()
    result = runner.invoke(app, [str(tmp_path), "--log-format", "invalid"])
    assert result.exit_code != 0
