"""Tests for the standalone profile validator CLI."""

import json
from pathlib import Path

import typer
from typer.testing import CliRunner

from codfreq.validate_profile import main


cli = typer.Typer()
cli.command()(main)


def _write_profile(path: Path, profile: dict) -> None:
    """Write *profile* as JSON to *path*."""

    path.write_text(json.dumps(profile), encoding="utf-8")


def test_validate_profile_cli_valid(tmp_path: Path) -> None:
    """CLI exits with code 0 for a valid profile."""

    profile = {
        "version": "20221213",
        "fragmentConfig": [{"fragmentName": "ref", "refSequence": "acgt"}],
        "sequenceAssemblyConfig": [],
    }
    prof = tmp_path / "prof.json"
    _write_profile(prof, profile)

    runner = CliRunner()
    result = runner.invoke(cli, [str(prof)])  # type: ignore[arg-type]
    assert result.exit_code == 0
    assert "Profile is valid" in result.stdout


def test_validate_profile_cli_invalid(tmp_path: Path) -> None:
    """CLI reports errors for an invalid profile."""

    profile = {
        "version": "20221213",
        "fragmentConfig": [],
        "sequenceAssemblyConfig": [{"trim": [1]}],
    }
    prof = tmp_path / "prof.json"
    _write_profile(prof, profile)

    runner = CliRunner()
    result = runner.invoke(cli, [str(prof)])  # type: ignore[arg-type]
    assert result.exit_code != 0
    assert "trim" in result.stdout
