"""Helpers for patching the optional PostAlign dependency."""

from __future__ import annotations

import sys
from contextlib import contextmanager
from types import ModuleType
from typing import Any, Iterator
from unittest.mock import patch


@contextmanager
def mock_postalign() -> Iterator[None]:
    """Patch ``postalign`` modules with lightweight stand-ins.

    These mocks provide just enough functionality for the parts of CodFreq
    exercised in tests while avoiding a runtime dependency on the actual
    PostAlign package.
    """
    postalign = ModuleType("postalign")
    postalign.__path__ = []  # type: ignore[attr-defined]
    utils = ModuleType("postalign.utils")
    utils.__path__ = []  # type: ignore[attr-defined]
    processors = ModuleType("postalign.processors")
    processors.__path__ = []  # type: ignore[attr-defined]
    models = ModuleType("postalign.models")
    models.__path__ = []  # type: ignore[attr-defined]

    cython = ModuleType("cython")

    def _decorator(*dargs: Any, **dkwargs: Any) -> Any:  # pragma: no cover
        if dargs and callable(dargs[0]):
            return dargs[0]
        return lambda func: func

    cython.ccall = _decorator  # type: ignore[attr-defined]
    cython.inline = _decorator  # type: ignore[attr-defined]
    cython.returns = lambda *a, **k: _decorator  # type: ignore[attr-defined]
    cython.void = None  # type: ignore[attr-defined]
    cython.cfunc = _decorator  # type: ignore[attr-defined]

    pysam = ModuleType("pysam")

    class AlignedSegment:  # pragma: no cover - placeholder
        """Minimal ``pysam.AlignedSegment`` stand-in."""

        def __init__(self, query_name: str | None = None) -> None:
            self.query_name = query_name

    class AlignmentFile:  # pragma: no cover - placeholder
        """Stub for ``pysam.AlignmentFile`` with a ``mapped`` count."""

        mapped = 0

        # pragma: no cover - stub
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            return

        def __enter__(self) -> "AlignmentFile":  # pragma: no cover
            return self

        def __exit__(self, *exc: Any) -> None:  # pragma: no cover
            return None

    pysam.AlignedSegment = AlignedSegment  # type: ignore[attr-defined]
    pysam.AlignmentFile = AlignmentFile  # type: ignore[attr-defined]

    def group_by_codons(
        seq1: bytearray, seq2: bytearray
    ) -> tuple[list[bytearray], list[bytearray]]:
        """Split two sequences into codon-sized chunks."""

        def chunk(seq: bytearray) -> list[bytearray]:
            return [seq[i:i + 3] for i in range(0, len(seq), 3)]

        return chunk(seq1), chunk(seq2)

    utils_group = ModuleType("postalign.utils.group_by_codons")
    utils_group.group_by_codons = group_by_codons  # type: ignore[attr-defined]

    def codon_align(
        refseq: Any, queryseq: Any, **_kwargs: Any
    ) -> tuple[Any, Any]:
        """Mutate ``queryseq`` to emulate codon realignment."""
        if hasattr(queryseq, "seqtext") and len(queryseq.seqtext) >= 3:
            queryseq.seqtext[:3] = b"TTT"
        return refseq, queryseq

    def parse_gap_placement_score(_text: str) -> dict:
        """Return an empty gap-placement score map."""
        return {}

    processors_codon = ModuleType("postalign.processors.codon_alignment")
    processors_codon.codon_align = codon_align  # type: ignore[attr-defined]
    processors_codon.parse_gap_placement_score = (
        parse_gap_placement_score  # type: ignore[attr-defined]
    )

    class NAPosition:
        """Minimal representation of a nucleotide position."""

        @staticmethod
        def init_from_bytes(b: bytearray) -> bytearray:
            """Return a mutable copy of *b*."""
            return bytearray(b)

        @staticmethod
        def as_bytes(seq: bytearray) -> bytes:
            """Return ``seq`` as immutable bytes."""
            return bytes(seq)

    class Sequence:
        """Simplified sequence container mimicking PostAlign's API."""

        def __init__(
            self,
            header: str,
            description: str,
            seqtext: bytearray,
            seqid: int,
            seqtype: Any,
            abs_seqstart: int,
            skip_invalid: bool,
        ) -> None:
            self.header = header
            self.description = description
            self.seqtext = seqtext
            self.seqid = seqid
            self.seqtype = seqtype
            self.abs_seqstart = abs_seqstart
            self.skip_invalid = skip_invalid

    models_seq = ModuleType("postalign.models.sequence")
    models_seq.NAPosition = NAPosition  # type: ignore[attr-defined]
    models_seq.Sequence = Sequence  # type: ignore[attr-defined]

    modules = {
        "postalign": postalign,
        "postalign.utils": utils,
        "postalign.utils.group_by_codons": utils_group,
        "postalign.processors": processors,
        "postalign.processors.codon_alignment": processors_codon,
        "postalign.models": models,
        "postalign.models.sequence": models_seq,
        "cython": cython,
        "pysam": pysam,
    }

    with patch.dict(sys.modules, modules):
        yield
