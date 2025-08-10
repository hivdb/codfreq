"""Tests for :mod:`codfreq.sam2consensus`."""

from __future__ import annotations

import json
import sys
import types
from typing import Any
from unittest.mock import MagicMock, mock_open, patch

# Stub cython decorators
cython_stub = types.ModuleType("cython")


def _decorator(*dargs: Any, **dkwargs: Any) -> Any:  # pragma: no cover
    if dargs and callable(dargs[0]):
        return dargs[0]
    return lambda func: func


cython_stub.cfunc = _decorator  # type: ignore[attr-defined]
cython_stub.ccall = _decorator  # type: ignore[attr-defined]
cython_stub.inline = _decorator  # type: ignore[attr-defined]
cython_stub.returns = lambda *a, **k: _decorator  # type: ignore[attr-defined]
sys.modules["cython"] = cython_stub

from codfreq.sam2consensus import (  # noqa: E402
    create_untrans_region_consensus,
    make_consensus,
    sam2consensus,
)
from codfreq.codfreq_types import NARegionConfig, Profile  # noqa: E402


def test_make_consensus_handles_gaps_and_insertions() -> None:
    """Missing bases yield gaps and insertions are appended."""

    region: NARegionConfig = {
        "name": "gag",
        "fromFragment": "frag",
        "refStart": 1,
        "refEnd": 2,
    }
    nacons_lookup: dict[tuple[int, int], int] = {
        (2, 0): ord("C"),
        (2, 1): ord("T"),
    }
    result = make_consensus(nacons_lookup, region)
    assert result == {
        "name": "gag",
        "refStart": 1,
        "refEnd": 2,
        "consensus": "NCT",
    }


@patch("codfreq.sam2consensus.get_posnas_in_genome_region")
def test_sam2consensus_keeps_frequent_insertions(
    mock_get_posnas: MagicMock,
) -> None:
    """Insertions are kept only when more than half of reads contain them."""

    mock_get_posnas.return_value = [
        (
            None,
            [
                (1, 0, ord("A"), 0),
                (1, 1, ord("T"), 0),
                (2, 0, ord("C"), 0),
            ],
        ),
        (
            None,
            [
                (1, 0, ord("A"), 0),
                (1, 1, ord("T"), 0),
                (2, 0, ord("C"), 0),
                (2, 1, ord("T"), 0),
            ],
        ),
        (
            None,
            [
                (1, 0, ord("A"), 0),
                (2, 0, ord("C"), 0),
            ],
        ),
    ]
    region: NARegionConfig = {
        "name": "gag",
        "fromFragment": "frag",
        "refStart": 1,
        "refEnd": 2,
    }
    result = sam2consensus("sample.sam", region)
    assert result == {
        "name": "gag",
        "refStart": 1,
        "refEnd": 2,
        "consensus": "ATC",
    }
    mock_get_posnas.assert_called_once_with(
        "sample.sam", ref_name="frag", ref_start=1, ref_end=2
    )


@patch("codfreq.sam2consensus.sam2consensus")
@patch("codfreq.sam2consensus.name_bamfile")
def test_create_untrans_region_consensus_writes_results(
    mock_name_bamfile: MagicMock,
    mock_sam2consensus: MagicMock,
) -> None:
    """Fragments without parents produce consensus JSON."""

    profile: Profile = {
        "version": "1",
        "fragmentConfig": [
            {"fragmentName": "F1"},
            {"fragmentName": "F2", "fromFragment": "other"},
        ],
        "sequenceAssemblyConfig": [
            {
                "name": "R1",
                "geneName": None,
                "fromFragment": "F1",
                "refStart": 1,
                "refEnd": 2,
            },
            {
                "name": "R2",
                "geneName": None,
                "fromFragment": "F2",
                "refStart": 1,
                "refEnd": 2,
            },
            {
                "name": None,
                "geneName": None,
                "fromFragment": "F1",
                "refStart": 1,
                "refEnd": 2,
            },
        ],
    }
    mock_name_bamfile.return_value = "file.bam"
    result_cons = {
        "name": "R1",
        "refStart": 1,
        "refEnd": 2,
        "consensus": "AT",
    }
    mock_sam2consensus.return_value = result_cons
    m = mock_open()
    with patch("codfreq.sam2consensus.open", m, create=True):
        create_untrans_region_consensus("seq1", profile)
    mock_name_bamfile.assert_called_once_with("seq1", "F1", is_trimmed=True)
    mock_sam2consensus.assert_called_once_with(
        "file.bam",
        {"name": "R1", "fromFragment": "F1", "refStart": 1, "refEnd": 2},
    )
    m.assert_called_once_with("seq1.untrans.json", "w")
    written = "".join(call.args[0] for call in m().write.call_args_list)
    assert json.loads(written) == [result_cons]
