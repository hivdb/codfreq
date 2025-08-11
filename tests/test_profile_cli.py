"""Tests for the profile CLI."""

import json
from pathlib import Path

import questionary  # type: ignore[import-not-found]
from typer.testing import CliRunner
from unittest.mock import patch

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
        "sequenceAssemblyConfig": [{"geneName": 1}],
    }
    prof = tmp_path / "prof.json"
    _write_profile(prof, profile)
    runner = CliRunner()
    result = runner.invoke(app, ["validate", str(prof)])
    assert result.exit_code != 0
    assert "geneName" in result.stdout


def test_create_profile(tmp_path: Path) -> None:
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

    with (
        patch.object(
            questionary, "text", side_effect=lambda m, **_: Prompt(m)
        ),
        patch.object(
            questionary, "confirm", side_effect=lambda m: Prompt(m)
        ),
        patch.object(
            questionary, "checkbox", side_effect=lambda *a, **k: Prompt("cb")
        ),
        patch.object(
            questionary, "select", side_effect=lambda *a, **k: Prompt("sel")
        ),
    ):
        runner = CliRunner()
        out = tmp_path / "new.json"
        result = runner.invoke(app, ["create", str(out)])
        assert result.exit_code == 0

    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["fragmentConfig"][0]["fragmentName"] == "ref"
    assert data["sequenceAssemblyConfig"][0]["refStart"] == 1
    Profile.model_validate(data)


def test_create_profile_genbank(tmp_path: Path) -> None:
    """Create a profile from a GenBank accession with gene suggestions."""

    record = profile_module.GenBankRecord(
        accession="TEST",
        sequence="ACGTACGTACGTACGTACGT",
        features=[
            profile_module.GeneFeature(name="gene1", ranges=[(1, 10)]),
            profile_module.GeneFeature(
                name="gene2", ranges=[(11, 15), (16, 20)]
            ),
        ],
    )

    answers = iter([
        "20221213",
        "TEST",
        "ref",
        "",
        True,
        True,
    ])

    class Prompt:
        def __init__(self, _msg: str) -> None:
            self._msg = _msg

        def ask(self) -> object:  # pragma: no cover - trivial
            return next(answers)

    def fake_checkbox(*_a: object, **_k: object) -> object:
        class P:
            def ask(self) -> object:  # pragma: no cover - trivial
                return record.features

        return P()

    with (
        patch.object(profile_module, "_fetch_record", return_value=record),
        patch.object(
            questionary, "text", side_effect=lambda m, **_: Prompt(m)
        ),
        patch.object(
            questionary, "confirm", side_effect=lambda m: Prompt(m)
        ),
        patch.object(questionary, "checkbox", side_effect=fake_checkbox),
    ):
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


def test_auto_assembly_overlap_and_gap() -> None:
    """Overlapping fragments yield left/right trimmed assemblies."""

    fragments: list[dict[str, object]] = [
        {"fragmentName": "ref", "refSequence": "A" * 25},
        {
            "fragmentName": "g1",
            "fromFragment": "ref",
            "geneName": "g1",
            "refRanges": [(1, 5)],
        },
        {
            "fragmentName": "g2",
            "fromFragment": "ref",
            "geneName": "g2",
            "refRanges": [(8, 15)],
        },
        {
            "fragmentName": "g3",
            "fromFragment": "ref",
            "geneName": "g3",
            "refRanges": [(14, 20)],
        },
    ]
    options = profile_module._auto_assembly_options(fragments)
    assert len(options) == 2
    left, right = options
    assert left[3]["trim"] == [[1, 2]]
    assert right[2]["trim"] == [[7, 8]]


def test_auto_assembly_no_main() -> None:
    """No main fragment produces no assembly options."""

    frags: list[dict[str, object]] = [
        {
            "fragmentName": "g1",
            "geneName": "g1",
            "fromFragment": "ref",
            "refRanges": [(1, 5)],
        }
    ]
    assert profile_module._auto_assembly_options(frags) == []
