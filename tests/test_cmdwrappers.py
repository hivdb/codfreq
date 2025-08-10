from pathlib import Path
from typing import Iterator

import pytest

import codfreq.compress_codfreq as cc
from codfreq.cmdwrappers.base import execute
from codfreq.enums import LogFormat


def test_find_codfreq_untrans_pairs(tmp_path: Path) -> None:
    cf = tmp_path / "sample.codfreq"
    cf.write_text("content", encoding="utf-8")
    ut = tmp_path / "sample.untrans.json"
    ut.write_text("[]", encoding="utf-8")
    hidden = tmp_path / ".hidden.codfreq"
    hidden.write_text("ignore", encoding="utf-8")
    other = tmp_path / "random.txt"
    other.write_text("ignore", encoding="utf-8")
    pairs = cc.find_codfreq_untrans_pairs(tmp_path)
    assert pairs == [(str(cf), str(ut))]


def test_find_codfreq_untrans_pairs_untrans_first(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    filenames = ["sample.untrans.json", "sample.codfreq"]

    def fake_walk(
        _path: Path, followlinks: bool = True
    ) -> Iterator[tuple[str, list[str], list[str]]]:
        yield (str(tmp_path), [], filenames)

    monkeypatch.setattr(cc.os, "walk", fake_walk)
    pairs = cc.find_codfreq_untrans_pairs(tmp_path)
    cf = tmp_path / "sample.codfreq"
    ut = tmp_path / "sample.untrans.json"
    assert pairs == [(str(cf), str(ut))]


def test_find_codfreq_untrans_pairs_codfreq_first(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """CodFreq files preceding untranslated files are paired."""

    filenames = ["sample.codfreq", "sample.untrans.json"]

    def fake_walk(
        _path: Path, followlinks: bool = True
    ) -> Iterator[tuple[str, list[str], list[str]]]:
        yield (str(tmp_path), [], filenames)

    monkeypatch.setattr(cc.os, "walk", fake_walk)
    pairs = cc.find_codfreq_untrans_pairs(tmp_path)
    cf = tmp_path / "sample.codfreq"
    ut = tmp_path / "sample.untrans.json"
    assert pairs == [(str(cf), str(ut))]


def test_compress_codfreq_logs_and_writes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    cf = tmp_path / "sample.codfreq"
    cf.write_text("header\n", encoding="utf-8")
    ut = tmp_path / "sample.untrans.json"
    ut.write_text(
        '[{"name": "UTR", "refStart": 1, "refEnd": 2, "consensus": "AT"}]',
        encoding="utf-8",
    )

    def fake_compress(
        data: bytes, compresslevel: int = 9, *, mtime: float | None = None
    ) -> bytes:
        return b"z" + bytes(data)

    monkeypatch.setattr(cc.pigz, "compress", fake_compress)
    printed: list[str] = []
    monkeypatch.setattr(cc.rich, "print", lambda msg: printed.append(msg))

    cc.compress_codfreq(tmp_path, log_format=LogFormat.json)
    out = capsys.readouterr().out
    assert '"op": "compress-codfreq"' in out
    gz = tmp_path / "sample.codfreq.gz"
    assert gz.read_bytes().startswith(b"z")
    gz.unlink()

    cc.compress_codfreq(tmp_path)
    assert printed == [f"Create {str(cf)}.gz"]


def test_execute_success_and_failure() -> None:
    import typer
    out, err = execute(["bash", "-c", "echo hello"])
    assert out.strip() == "hello"
    with pytest.raises(typer.Abort):
        execute(["bash", "-c", "exit 1"])
