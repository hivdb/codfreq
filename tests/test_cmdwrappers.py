from pathlib import Path
import sys

import pytest

# Use shimmed typer for tests when real package is unavailable
shims_path = Path(__file__).resolve().parent / "shims"
sys.path.insert(0, str(shims_path))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from codfreq.compress_codfreq import find_codfreq_untrans_pairs
from codfreq.cmdwrappers.base import execute


def test_find_codfreq_untrans_pairs(tmp_path: Path) -> None:
    cf = tmp_path / "sample.codfreq"
    cf.write_text("content", encoding="utf-8")
    ut = tmp_path / "sample.untrans.json"
    ut.write_text("[]", encoding="utf-8")
    pairs = find_codfreq_untrans_pairs(tmp_path)
    assert pairs == [(str(cf), str(ut))]


def test_execute_success_and_failure() -> None:
    import typer
    out, err = execute(["bash", "-c", "echo hello"])
    assert out.strip() == "hello"
    with pytest.raises(typer.Abort):
        execute(["bash", "-c", "exit 1"])
