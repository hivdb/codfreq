from pathlib import Path
import gzip
from typing import Any

from typer.testing import CliRunner

from codfreq.compress_codfreq import app
from codfreq.cmdwrappers import pigz


def test_compress_codfreq_cli(tmp_path: Path, monkeypatch: Any) -> None:
    cf = tmp_path / "sample.codfreq"
    cf.write_text("content", encoding="utf-8")
    ut = tmp_path / "sample.untrans.json"
    ut.write_text("[]", encoding="utf-8")
    monkeypatch.setattr(pigz, "compress", lambda data: gzip.compress(data))
    runner = CliRunner()
    result = runner.invoke(app, [str(tmp_path)])
    assert result.exit_code == 0
    assert (tmp_path / "sample.codfreq.gz").exists()


def test_compress_cli_json_logging(tmp_path: Path, monkeypatch: Any) -> None:
    """The CLI logs JSON when requested."""

    cf = tmp_path / "sample.codfreq"
    cf.write_text("content", encoding="utf-8")
    ut = tmp_path / "sample.untrans.json"
    ut.write_text("[]", encoding="utf-8")
    monkeypatch.setattr(pigz, "compress", lambda data: gzip.compress(data))
    runner = CliRunner()
    result = runner.invoke(
        app, [str(tmp_path), "--log-format", "json"])
    assert result.exit_code == 0
    assert "compress-codfreq" in result.stdout
    assert (tmp_path / "sample.codfreq.gz").exists()


def test_compress_codfreq_cli_invalid_log_format(tmp_path: Path) -> None:
    runner = CliRunner()
    result = runner.invoke(app, [str(tmp_path), "--log-format", "invalid"])
    assert result.exit_code != 0
