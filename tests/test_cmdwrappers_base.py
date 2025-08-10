"""Tests for cmdwrapper base utilities."""

from unittest.mock import MagicMock, patch

import pytest
import typer

from codfreq.cmdwrappers import base


def test_execute_invokes_subprocess() -> None:
    """Commands run via ``execute`` capture stdout and stderr."""
    proc = MagicMock()
    proc.communicate.return_value = ("out", "err")
    proc.returncode = 0
    with (
        patch(
            "codfreq.cmdwrappers.base.Popen",
            return_value=proc,
        ) as popen_mock,
        patch("codfreq.cmdwrappers.base.raise_on_proc_error") as raise_mock,
    ):
        out, err = base.execute(["echo", "hi"])
    assert (out, err) == ("out", "err")
    popen_mock.assert_called_once_with(
        ["echo", "hi"],
        stdout=base.PIPE,
        stderr=base.PIPE,
        encoding="U8",
    )
    raise_mock.assert_called_once_with(proc, "err")


def test_raise_on_proc_error_raises() -> None:
    """Non-zero return codes trigger a ``typer.Abort``."""
    proc = MagicMock(returncode=1)
    with patch("rich.print") as rich_print, pytest.raises(typer.Abort):
        base.raise_on_proc_error(proc, "bad")
    rich_print.assert_called_once()


def test_raise_on_proc_error_noop() -> None:
    """Zero return codes result in no output or exception."""
    proc = MagicMock(returncode=0)
    with patch("rich.print") as rich_print:
        base.raise_on_proc_error(proc, "")
    rich_print.assert_not_called()


def test_decorator_registration() -> None:
    """Decorator helpers register functions for later lookup."""
    def func() -> None:
        return None

    wrapped_ref = base.refinit_func("demo")(func)
    assert base.get_refinit("demo") is wrapped_ref

    wrapped_align = base.align_func("demo2")(func)
    assert "demo2" in base.get_programs()
    assert base.get_align("demo2") is wrapped_align
