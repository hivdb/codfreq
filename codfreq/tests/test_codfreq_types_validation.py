"""Validation tests for configuration models."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from codfreq.codfreq_types import (
    CodonAlignmentConfig,
    DerivedFragmentConfig,
    GeneAssemblyConfig,
    RegionAssemblyConfig,
    Profile,
)


def test_codon_alignment_config_rejects_extra() -> None:
    """``CodonAlignmentConfig`` raises on unknown fields."""
    with pytest.raises(ValidationError):
        CodonAlignmentConfig(  # type: ignore[call-arg]
            relRefStart=1,
            relRefEnd=2,
            unknown=3,
        )


def test_derived_fragment_config_rejects_extra() -> None:
    """``DerivedFragmentConfig`` raises on unknown fields."""
    with pytest.raises(ValidationError):
        DerivedFragmentConfig(  # type: ignore[call-arg]
            fragmentName="frag",
            fromFragment="ref",
            refRanges=[(1, 2)],
            unknown=5,
        )


def test_gene_assembly_config_rejects_extra() -> None:
    """``GeneAssemblyConfig`` raises on unknown fields."""
    with pytest.raises(ValidationError):
        GeneAssemblyConfig(  # type: ignore[call-arg]
            geneName="a", unknown=1
        )


def test_region_assembly_config_rejects_extra() -> None:
    """``RegionAssemblyConfig`` raises on unknown fields."""
    with pytest.raises(ValidationError):
        RegionAssemblyConfig(  # type: ignore[call-arg]
            name="a", fromFragment="ref", refStart=1, refEnd=2, unknown=1
        )


def test_gene_assembly_trim_normalization() -> None:
    """Single integers become ``(n, n)`` ranges."""

    cfg = GeneAssemblyConfig(geneName="g", trim=[1, (5, 6)])
    assert cfg.trim == [(1, 1), (5, 6)]


def test_profile_rejects_extra() -> None:
    """``Profile`` raises on unknown fields."""
    with pytest.raises(ValidationError):
        Profile(  # type: ignore[call-arg]
            version="1",
            fragmentConfig=[],
            sequenceAssemblyConfig=[],
            unknown=1,
        )


def test_profile_assembly_contiguity() -> None:
    """Profile validator enforces contiguous assembly coverage."""

    main = {"fragmentName": "ref", "refSequence": "acgt"}
    valid = {
        "version": "1",
        "fragmentConfig": [main],
        "sequenceAssemblyConfig": [
            {"name": "ref", "fromFragment": "ref", "refStart": 1, "refEnd": 4}
        ],
    }
    Profile.model_validate(valid)
    invalid = {
        "version": "1",
        "fragmentConfig": [main],
        "sequenceAssemblyConfig": [
            {"name": "ref", "fromFragment": "ref", "refStart": 2, "refEnd": 4}
        ],
    }
    with pytest.raises(ValidationError):
        Profile.model_validate(invalid)
    truncated = {
        "version": "1",
        "fragmentConfig": [main],
        "sequenceAssemblyConfig": [
            {"name": "ref", "fromFragment": "ref", "refStart": 1, "refEnd": 3}
        ],
    }
    with pytest.raises(ValidationError):
        Profile.model_validate(truncated)


def test_profile_unknown_gene() -> None:
    """Unknown genes in assemblies raise an error."""

    main = {"fragmentName": "ref", "refSequence": "acgt"}
    profile = {
        "version": "1",
        "fragmentConfig": [main],
        "sequenceAssemblyConfig": [
            {"geneName": "g1"},
        ],
    }
    with pytest.raises(ValidationError) as err:
        Profile.model_validate(profile)
    assert "Unknown gene" in str(err.value)
