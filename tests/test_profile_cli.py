"""Tests for the profile CLI."""

import json
import io
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import cast

try:  # pragma: no cover - dependency shim
    import questionary  # type: ignore[import-not-found]
except ModuleNotFoundError:  # pragma: no cover - dependency shim
    questionary = SimpleNamespace(  # type: ignore[assignment]
        text=lambda *a, **k: None,
        confirm=lambda *a, **k: None,
        checkbox=lambda *a, **k: None,
        Choice=lambda *a, **k: SimpleNamespace(**k),
    )
    sys.modules["questionary"] = cast(ModuleType, questionary)

try:  # pragma: no cover - dependency shim
    from Bio import Entrez, SeqIO  # type: ignore[import-not-found]
except ModuleNotFoundError:  # pragma: no cover - dependency shim
    Entrez = SimpleNamespace(efetch=lambda *a, **k: None)
    SeqIO = SimpleNamespace(read=lambda *a, **k: None)
    bio_mod = SimpleNamespace(Entrez=Entrez, SeqIO=SeqIO)
    sys.modules["Bio"] = cast(ModuleType, bio_mod)
    sys.modules["Bio.Entrez"] = Entrez  # type: ignore[assignment]
    sys.modules["Bio.SeqIO"] = SeqIO  # type: ignore[assignment]
    sys.modules["Bio.SeqFeature"] = SimpleNamespace(
        CompoundLocation=type("CompoundLocation", (), {})
    )  # type: ignore[assignment]
    from Bio import Entrez, SeqIO  # type: ignore[import-not-found]

from typer.testing import CliRunner
from unittest.mock import patch

import codfreq.profile as profile_module
from codfreq.profile import app, validate_app
from codfreq.codfreq_types import (
    Profile,
    FragmentConfig,
    MainFragmentConfig,
    DerivedFragmentConfig,
    GeneAssemblyConfig,
    RegionAssemblyConfig,
)


def _write_profile(path: Path, profile: dict) -> None:
    """Write profile JSON to disk."""
    path.write_text(json.dumps(profile), encoding="utf-8")


def test_fetch_record_parses_features() -> None:
    """``_fetch_record`` extracts sequences and gene ranges."""

    feature_skip = SimpleNamespace(
        type="misc_feature",
        qualifiers={},
        location=SimpleNamespace(start=0, end=1),
    )
    feature_missing = SimpleNamespace(
        type="gene",
        qualifiers={},
        location=SimpleNamespace(start=1, end=2),
    )
    feature = SimpleNamespace(
        type="gene",
        qualifiers={"gene": ["g1"]},
        location=SimpleNamespace(start=2, end=4),
    )
    record = SimpleNamespace(
        seq="ACGT",
        features=[feature_skip, feature_missing, feature],
    )

    class Handle(io.StringIO):
        def __enter__(self) -> "Handle":
            return self

        def __exit__(self, *exc: object) -> None:
            self.close()

    with (
        patch.object(profile_module.Entrez, "efetch", return_value=Handle()),
        patch.object(profile_module.SeqIO, "read", return_value=record),
    ):
        result = profile_module._fetch_record("ACC", "a@b")
    assert result.accession == "ACC"
    assert result.sequence == "ACGT"
    assert result.features[0].name == "g1"
    assert result.features[0].ranges == [(3, 4)]


def test_validate_profile_valid(tmp_path: Path) -> None:
    """Standalone validator exits with code 0 for valid input."""
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
        "sequenceAssemblyConfig": [
            {
                "name": "ref",
                "fromFragment": "ref",
                "refStart": 1,
                "refEnd": 4,
            }
        ],
    }
    prof = tmp_path / "prof.json"
    _write_profile(prof, profile)
    runner = CliRunner()
    result = runner.invoke(validate_app, [str(prof)])
    assert result.exit_code == 0


def test_validate_profile_invalid(tmp_path: Path) -> None:
    """Standalone validator reports errors for invalid input."""
    profile = {
        "version": "20221213",
        "fragmentConfig": [{"fragmentName": "ref", "refSequence": "acgt"}],
        "sequenceAssemblyConfig": [{"geneName": 1}],
    }
    prof = tmp_path / "prof.json"
    _write_profile(prof, profile)
    runner = CliRunner()
    result = runner.invoke(validate_app, [str(prof)])
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
        True,
    ])

    class Prompt:
        def __init__(self, _msg: str) -> None:
            self._msg = _msg

        def ask(self) -> object:
            return next(answers)

    with (
        patch.object(
            questionary, "text", side_effect=lambda m, **_: Prompt(m)
        ),
        patch.object(
            questionary, "confirm", side_effect=lambda m: Prompt(m)
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
        "user@example.com",
        "ref",
        "",
        True,
    ])

    class Prompt:
        def __init__(self, _msg: str) -> None:
            self._msg = _msg

        def ask(self) -> object:
            return next(answers)

    def fake_checkbox(*_a: object, **_k: object) -> object:
        class P:
            def ask(self) -> object:
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


def test_prompt_genbank_fragments_no_features() -> None:
    """Records without features yield no fragments."""

    record = profile_module.GenBankRecord(
        accession="A", sequence="ACGT", features=[]
    )
    assert (
        profile_module._prompt_genbank_fragments(record, "ref") == []
    )


def test_prompt_manual_fragments() -> None:
    """Manual fragment prompts default gene name to fragment name."""

    answers = iter([
        "ref",
        "ACGT",
        "g1",
        "",
        "1-3,5,",
        "",
    ])

    class Prompt:
        def __init__(self, _m: str) -> None:
            self._m = _m

        def ask(self) -> object:
            return next(answers)

    with patch.object(
        questionary, "text", side_effect=lambda m, **_: Prompt(m)
    ):
        main, frags = profile_module._prompt_manual_fragments(None)
    assert main == "ref"
    assert isinstance(frags[0], MainFragmentConfig)
    derived = frags[1]
    assert isinstance(derived, DerivedFragmentConfig)
    assert derived.geneName == "g1"
    assert derived.refRanges == [(1, 3), (5, 5)]


def test_prompt_manual_assemblies() -> None:
    """Manual assembly prompts support gene and inter-gene regions."""

    frags: list[FragmentConfig] = [
        MainFragmentConfig(fragmentName="ref", refSequence="AAAAA"),
        DerivedFragmentConfig(
            fragmentName="g1",
            fromFragment="ref",
            geneName="g1",
            refRanges=[(1, 3)],
        ),
    ]

    answers = iter([
        True,  # add region
        True,  # is gene
        "g1",
        "1,3-4,",  # trim
        True,  # add region
        False,  # not gene
        "tail",
        "ref",
        "4",
        "5",
        False,  # stop
    ])

    class Prompt:
        def __init__(self, _m: str) -> None:
            self._m = _m

        def ask(self) -> object:
            return next(answers)

    with (
        patch.object(questionary, "confirm", side_effect=lambda m: Prompt(m)),
        patch.object(
            questionary, "text", side_effect=lambda m, **_: Prompt(m)
        ),
    ):
        assemblies = profile_module._prompt_manual_assemblies(frags)
    assert isinstance(assemblies[0], GeneAssemblyConfig)
    assert isinstance(assemblies[1], RegionAssemblyConfig)


def test_create_profile_rejects_suggestion(tmp_path: Path) -> None:
    """User can override suggested assembly with manual entries."""

    answers = iter([
        "20221213",
        "",  # no GenBank
        "ref",
        "ac",
        "g1",
        "g1",
        "1-2",
        "",  # stop fragments
        False,  # reject auto assembly
        True,  # add region
        True,  # gene
        "g1",
        "",  # trim
        False,  # stop
    ])

    class Prompt:
        def __init__(self, _m: str) -> None:
            self._m = _m

        def ask(self) -> object:
            return next(answers)

    with (
        patch.object(
            questionary, "text", side_effect=lambda m, **_: Prompt(m)
        ),
        patch.object(questionary, "confirm", side_effect=lambda m: Prompt(m)),
    ):
        runner = CliRunner()
        out = tmp_path / "override.json"
        result = runner.invoke(app, ["create", str(out)])
        assert result.exit_code == 0


def test_create_profile_manual_assembly_when_missing(tmp_path: Path) -> None:
    """Manual assembly is requested when no suggestion is available."""

    answers = iter([
        "20221213",
        "",  # no GenBank
        "ref",
        "acgt",
        "",  # stop fragments
        # no confirm because assemblies is empty
        True,
        False,
    ])

    class Prompt:
        def __init__(self, _m: str) -> None:
            self._m = _m

        def ask(self) -> object:
            return next(answers)

    with (
        patch.object(profile_module, "_auto_assembly", return_value=[]),
        patch.object(
            profile_module,
            "_prompt_manual_assemblies",
            return_value=[
                RegionAssemblyConfig(
                    name="ref", fromFragment="ref", refStart=1, refEnd=4
                )
            ],
        ),
        patch.object(
            questionary, "text", side_effect=lambda m, **_: Prompt(m)
        ),
        patch.object(questionary, "confirm", side_effect=lambda m: Prompt(m)),
    ):
        runner = CliRunner()
        out = tmp_path / "manual.json"
        result = runner.invoke(app, ["create", str(out)])
        assert result.exit_code == 0


def test_create_profile_validation_failure(tmp_path: Path) -> None:
    """Invalid assembly data triggers validation error."""

    answers = iter([
        "20221213",
        "",  # no GenBank
        "ref",
        "acgt",
        "",  # stop fragments
        True,
    ])

    class Prompt:
        def __init__(self, _m: str) -> None:
            self._m = _m

        def ask(self) -> object:
            return next(answers)

    bad_assembly = [
        RegionAssemblyConfig(
            name="ref", fromFragment="ref", refStart=2, refEnd=4
        )
    ]

    with (
        patch.object(
            profile_module, "_auto_assembly", return_value=bad_assembly
        ),
        patch.object(
            questionary, "text", side_effect=lambda m, **_: Prompt(m)
        ),
        patch.object(questionary, "confirm", side_effect=lambda m: Prompt(m)),
    ):
        runner = CliRunner()
        out = tmp_path / "bad.json"
        result = runner.invoke(app, ["create", str(out)])
        assert result.exit_code != 0


def test_auto_assembly_overlap_and_gap() -> None:
    """Overlapping fragments yield left-trimmed assemblies."""

    fragments: list[FragmentConfig] = [
        MainFragmentConfig(fragmentName="ref", refSequence="A" * 25),
        DerivedFragmentConfig(
            fragmentName="g1",
            fromFragment="ref",
            geneName="g1",
            refRanges=[(1, 5)],
        ),
        DerivedFragmentConfig(
            fragmentName="g2",
            fromFragment="ref",
            geneName="g2",
            refRanges=[(8, 15)],
        ),
        DerivedFragmentConfig(
            fragmentName="g3",
            fromFragment="ref",
            geneName="g3",
            refRanges=[(14, 20)],
        ),
    ]
    assemblies = profile_module._auto_assembly(fragments)
    assert isinstance(assemblies[3], GeneAssemblyConfig)
    assert assemblies[3].trim == [(1, 2)]


def test_auto_assembly_no_main() -> None:
    """No main fragment produces no assembly options."""

    frags: list[FragmentConfig] = [
        DerivedFragmentConfig(
            fragmentName="g1",
            fromFragment="ref",
            geneName="g1",
            refRanges=[(1, 5)],
        )
    ]
    assert profile_module._auto_assembly(frags) == []
