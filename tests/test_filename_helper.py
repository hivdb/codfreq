"""Tests for :mod:`codfreq.filename_helper`."""

from codfreq.filename_helper import (
    name_bamfile,
    name_codfreq,
    name_file,
    replace_ext,
    suggest_pair_name,
)


def test_replace_ext_basic() -> None:
    """``replace_ext`` swaps extensions and strips the old one."""

    assert replace_ext("sample.fastq", ".bam") == "sample.bam"
    renamed = replace_ext(
        "/path/sample.fastq", ".sam", fromext=".fastq", name_only=True
    )
    assert renamed == "sample.sam"


def test_name_helpers() -> None:
    """Check basic file naming helpers."""

    assert name_codfreq("run1") == "run1.codfreq"
    assert name_bamfile("run1", "ref") == "run1.ref.bam"
    assert name_bamfile("run1", "ref", is_trimmed=False) == "run1.ref.orig.bam"


def test_suggest_and_name_file() -> None:
    """``suggest_pair_name`` drives ``name_file``."""

    fnpair = ("reads_R1.fastq", "reads_R2.fastq")
    pattern = ("_", 1, 0, 0)
    assert suggest_pair_name(fnpair, pattern) == "reads"
    assert name_file(fnpair, pattern, ".bam") == "reads.bam"
