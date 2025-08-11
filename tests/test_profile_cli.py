"""Tests for the profile CLI."""

import json
from pathlib import Path

import questionary  # type: ignore[import-not-found]
from pytest import MonkeyPatch
from typer.testing import CliRunner

import codfreq.profile as profile_module
from codfreq.profile import app
from codfreq.codfreq_types import Profile


def _write_profile(path: Path, profile: dict) -> None:
    """Write profile JSON to disk."""
    path.write_text(json.dumps(profile), encoding="utf-8")


def test_validate_profile_valid(tmp_path: Path) -> None:
    """CLI exits with code 0 for a valid profile."""
    profile = {
        "version": "20221213",
        "fragmentConfig": [
            {"fragmentName": "ref", "refSequence": "acgt"},
            {
                "fragmentName": "frag1",
                "fromFragment": "ref",
                "refRanges": [[1, 4]],
            },
        ],
        "sequenceAssemblyConfig": [],
    }
    prof = tmp_path / "prof.json"
    _write_profile(prof, profile)
    runner = CliRunner()
    result = runner.invoke(app, ["validate", str(prof)])
    assert result.exit_code == 0
    assert "Profile is valid" in result.stdout


def test_validate_profile_invalid(tmp_path: Path) -> None:
    """CLI reports errors for an invalid profile."""
    profile = {
        "version": "20221213",
        "fragmentConfig": [{"fragmentName": "ref", "refSequence": "acgt"}],
        "sequenceAssemblyConfig": [{"trim": [1]}],
    }
    prof = tmp_path / "prof.json"
    _write_profile(prof, profile)
    runner = CliRunner()
    result = runner.invoke(app, ["validate", str(prof)])
    assert result.exit_code != 0
    assert "trim" in result.stdout


def test_create_profile(tmp_path: Path, monkeypatch: MonkeyPatch) -> None:
    """Interactively create a minimal profile."""

    answers = iter([
        "20221213",
        "",  # no GenBank accession
        "ref",
        "acgt",
        "",
        False,
        True,
        "region1",
        "ref",
        "1",
        "2",
        False,
    ])

    class Prompt:
        def __init__(self, _msg: str) -> None:
            self._msg = _msg

        def ask(self) -> object:  # pragma: no cover - trivial
            return next(answers)

    monkeypatch.setattr(questionary, "text", lambda msg, **_: Prompt(msg))
    monkeypatch.setattr(questionary, "confirm", lambda msg: Prompt(msg))
    monkeypatch.setattr(questionary, "checkbox", lambda *a, **k: Prompt("cb"))

    runner = CliRunner()
    out = tmp_path / "new.json"
    result = runner.invoke(app, ["create", str(out)])
    assert result.exit_code == 0

    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["fragmentConfig"][0]["fragmentName"] == "ref"
    assert data["sequenceAssemblyConfig"][0]["refStart"] == 1
    Profile.model_validate(data)


def test_create_profile_genbank(
    tmp_path: Path, monkeypatch: MonkeyPatch
) -> None:
    """Create a profile from a GenBank accession with gene suggestions."""

    gb_record = (
        "LOCUS       TEST        20 bp    DNA     linear   01-JAN-2000\n"
        "DEFINITION  Test sequence\n"
        "ACCESSION   TEST\n"
        "VERSION     TEST.1\n"
        "FEATURES             Location/Qualifiers\n"
        "     gene            1..10\n"
        "                     /gene=\"gene1\"\n"
        "     gene            join(11..15,16..20)\n"
        "                     /gene=\"gene2\"\n"
        "ORIGIN\n"
        "        1 acgtacgtacgtacgtacgt\n"
        "//\n"
    )

    monkeypatch.setattr(
        profile_module, "_download_genbank", lambda _acc: gb_record
    )
    record = profile_module._parse_genbank(gb_record, "TEST")

    answers = iter([
        "20221213",
        "TEST",
        "ref",
        "",
        True,
    ])

    class Prompt:
        def __init__(self, _msg: str) -> None:
            self._msg = _msg

        def ask(self) -> object:  # pragma: no cover - trivial
            return next(answers)

    monkeypatch.setattr(questionary, "text", lambda msg, **_: Prompt(msg))
    monkeypatch.setattr(questionary, "confirm", lambda msg: Prompt(msg))

    def fake_checkbox(*_a: object, **_k: object) -> object:
        class CB:
            def ask(self) -> object:  # pragma: no cover - trivial
                return record.features

        return CB()

    monkeypatch.setattr(questionary, "checkbox", fake_checkbox)

    runner = CliRunner()
    out = tmp_path / "gb.json"
    result = runner.invoke(app, ["create", str(out)])
    assert result.exit_code == 0

    data = json.loads(out.read_text(encoding="utf-8"))
    assert len(data["fragmentConfig"]) == 3
    gene2 = next(
        f for f in data["fragmentConfig"] if f["fragmentName"] == "gene2"
    )
    assert gene2["refRanges"] == [[11, 15], [16, 20]]
    Profile.model_validate(data)
