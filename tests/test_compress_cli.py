from pathlib import Path
import gzip

from typer.testing import CliRunner

from codfreq.compress_codfreq import app
from codfreq.cmdwrappers import pigz


def test_compress_codfreq_cli(tmp_path: Path, monkeypatch) -> None:
    cf = tmp_path / "sample.codfreq"
    cf.write_text("content", encoding="utf-8")
    ut = tmp_path / "sample.untrans.json"
    ut.write_text("[]", encoding="utf-8")
    monkeypatch.setattr(pigz, "compress", lambda data: gzip.compress(data))
    runner = CliRunner()
    result = runner.invoke(app, [str(tmp_path)])
    assert result.exit_code == 0
    assert (tmp_path / "sample.codfreq.gz").exists()
