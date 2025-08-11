"""Validation tests for configuration models."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from codfreq.codfreq_types import (
    CodonAlignmentConfig,
    DerivedFragmentConfig,
    SequenceAssemblyConfig,
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


def test_sequence_assembly_config_rejects_extra() -> None:
    """``SequenceAssemblyConfig`` raises on unknown fields."""
    with pytest.raises(ValidationError):
        SequenceAssemblyConfig(name="a", unknown=1)  # type: ignore[call-arg]


def test_sequence_assembly_trim_normalization() -> None:
    """Single integers become ``(n, n)`` ranges."""

    cfg = SequenceAssemblyConfig(trim=[1, (5, 6)])
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
