"""Tests for :mod:`codfreq.sam2codfreq`."""

from collections import Counter
from typing import cast

from unittest.mock import MagicMock, patch

import codfreq.sam2codfreq as s2c
from codfreq.codfreq_types import (
    Profile,
    MainFragmentConfig,
    DerivedFragmentConfig,
)


def test_build_fragment_intervals_and_get_ref_ranges() -> None:
    fragments = [
        {"refRanges": [(1, 3), (5, 7)], "fragmentName": "fragA"},
        {"refRanges": [(8, 9)], "fragmentName": "fragB"},
    ]
    assert s2c.build_fragment_intervals(fragments) == [
        ([(1, 3), (5, 7)], "fragA"),
        ([(8, 9)], "fragB"),
    ]

    assert s2c.get_ref_ranges({"refRanges": [(1, 3)]}) == [(1, 3)]
    assert s2c.get_ref_ranges({"refStart": 4, "refEnd": 6}) == [(4, 6)]


def test_get_ref_fragments() -> None:
    profile = Profile.model_validate(
        {
            "fragmentConfig": [
                {"fragmentName": "refA", "refSequence": "AAA"},
                {
                    "fragmentName": "fragA",
                    "fromFragment": "refA",
                    "geneName": "geneX",
                    "refRanges": [(1, 3)],
                    "codonAlignment": [
                        {
                            "relRefStart": 1,
                            "relRefEnd": 3,
                            "windowSize": 3,
                            "minGapDistance": 5,
                            "relGapPlacementScore": "0:0-1=1",
                        }
                    ],
                },
                {
                    "fragmentName": "fragB",
                    "fromFragment": "refA",
                    "geneName": "geneY",
                    "refRanges": [(4, 6)],
                    "codonAlignment": False,
                },
            ]
        }
    )
    refs, lookup = s2c.get_ref_fragments(profile)
    assert refs[0][2][0].codonAlignment[0].relRefEnd == 3
    assert refs[0][2][0].codonAlignment[0].minGapDistance == 5
    assert (
        refs[0][2][0].codonAlignment[0].relGapPlacementScore == "0:0-1=1"
    )
    assert refs[0][2][1].codonAlignment is False
    assert lookup["fragB"] == [("geneY", 0)]


def test_to_codon_counter_by_fragpos_and_get_codonfreq() -> None:
    codon_counter = Counter({("fragA", 1, "AAA"): 2, ("fragA", 1, "CCC"): 1})
    fragpos = s2c.to_codon_counter_by_fragpos(codon_counter)
    assert fragpos == {("fragA", 1): Counter({"AAA": 2, "CCC": 1})}

    qualities = {("fragA", 1): Counter({"AAA": 50, "CCC": 20})}
    lookup = {"fragA": [("geneX", 0)]}
    rows = s2c.get_codonfreq(fragpos, qualities, lookup)
    assert rows == [
        {
            "gene": "geneX",
            "position": 1,
            "total": 3,
            "codon": "AAA",
            "count": 2,
            "total_quality_score": 50,
        },
        {
            "gene": "geneX",
            "position": 1,
            "total": 3,
            "codon": "CCC",
            "count": 1,
            "total_quality_score": 20,
        },
    ]


def test_sam2codfreq_between() -> None:
    """Aggregate codons from positional reads."""

    iter_mock = MagicMock(
        return_value=iter(
            [
                (None, [("fragA", 1, "AAA", 30)]),
                (
                    None,
                    [("fragA", 1, "AAA", 20), ("fragA", 2, "CCC", 10)],
                ),
            ]
        )
    )

    with patch.object(s2c, "iter_poscodons", iter_mock):
        stat, qual, num_row = s2c.sam2codfreq_between(
            "file.bam", 0, 10, [([(1, 3)], "fragA")]
        )

    assert stat == Counter({("fragA", 1, "AAA"): 2, ("fragA", 2, "CCC"): 1})
    assert qual == Counter({("fragA", 1, "AAA"): 50, ("fragA", 2, "CCC"): 10})
    assert num_row == 2


def test_sam2codfreq_all() -> None:
    """Process all fragments and convert to rows."""

    with patch.object(
        s2c,
        "sam2codfreq",
        return_value=(
            {("fragA", 1): Counter({"AAA": 1})},
            {("fragA", 1): Counter({"AAA": 30})},
        ),
    ), patch(
        "codfreq.sam2codfreq.name_bamfile", return_value="file.bam"
    ):
        profile = Profile.model_validate(
            {
                "fragmentConfig": [
                    {"fragmentName": "refA", "refSequence": "AAA"},
                    {
                        "fragmentName": "fragA",
                        "fromFragment": "refA",
                        "geneName": "geneX",
                        "refRanges": [(1, 3)],
                    },
                ]
            }
        )

        rows = s2c.sam2codfreq_all("sample", (None, None), profile, workers=1)

    assert rows == [
        {
            "gene": "geneX",
            "position": 1,
            "total": 1,
            "codon": "AAA",
            "count": 1,
            "total_quality_score": 30,
        }
    ]


def test_sam2codfreq_accumulates_and_reports_progress() -> None:
    """Multiprocessing results are merged and progress is reported."""

    class DummyExecutor:
        """Context manager yielding preset ``map`` results."""

        def __init__(self, _workers: int) -> None:
            return

        def __enter__(self) -> "DummyExecutor":
            return self

        def __exit__(self, *exc: object) -> bool:
            return False

        def map(self, func, *args):  # type: ignore[no-untyped-def]
            return iter(
                [
                    (
                        Counter({("fragA", 1, "AAA"): 1}),
                        Counter({("fragA", 1, "AAA"): 30}),
                        1,
                    ),
                    (
                        Counter({("fragA", 2, "CCC"): 1}),
                        Counter({("fragA", 2, "CCC"): 20}),
                        1,
                    ),
                ]
            )

    class DummyAlignmentFile:
        mapped = 2

        def __enter__(self) -> "DummyAlignmentFile":
            return self

        def __exit__(self, *exc: object) -> None:
            return None

    progress = MagicMock()

    ref = MainFragmentConfig(fragmentName="refA", refSequence="AAA")
    fragments = [
        DerivedFragmentConfig(
            fragmentName="fragA",
            fromFragment="refA",
            refRanges=[(1, 3)],
        )
    ]

    with (
        patch(
            "codfreq.sam2codfreq.pysam.AlignmentFile",
            return_value=DummyAlignmentFile(),
        ),
        patch("codfreq.sam2codfreq.tqdm", return_value=progress),
        patch("codfreq.sam2codfreq.ProcessPoolExecutor", DummyExecutor),
        patch(
            "codfreq.sam2codfreq.chunked_samfile",
            return_value=[(0, 1), (1, 2)],
        ),
        patch(
            "codfreq.sam2codfreq.codonalign_consensus",
            side_effect=lambda *a: a[:2],
        ),
    ):
        codonstat, qualities = s2c.sam2codfreq(
            "file.bam", ref, cast(list, fragments), workers=1
        )

    assert codonstat == {
        ("fragA", 1): Counter({"AAA": 1}),
        ("fragA", 2): Counter({"CCC": 1}),
    }
    assert qualities == {
        ("fragA", 1): Counter({"AAA": 30}),
        ("fragA", 2): Counter({"CCC": 20}),
    }
    progress.set_description.assert_called_once()
    assert progress.update.call_count == 2
    progress.close.assert_called_once()


def test_sam2codfreq_supports_json_log_format() -> None:
    """JSON log format initializes :class:`JsonProgress`."""

    class DummyExecutor:
        def __init__(self, _workers: int) -> None:
            return

        def __enter__(self) -> "DummyExecutor":
            return self

        def __exit__(self, *exc: object) -> bool:
            return False

        def map(self, func, *args):  # type: ignore[no-untyped-def]
            return iter([(Counter(), Counter(), 0)])

    class DummyAlignmentFile:
        mapped = 0

        def __enter__(self) -> "DummyAlignmentFile":
            return self

        def __exit__(self, *exc: object) -> None:
            return None

    json_progress = MagicMock()

    ref = MainFragmentConfig(fragmentName="refA", refSequence="AAA")
    fragments = [
        DerivedFragmentConfig(
            fragmentName="fragA",
            fromFragment="refA",
            refRanges=[(1, 3)],
        )
    ]

    with (
        patch(
            "codfreq.sam2codfreq.pysam.AlignmentFile",
            return_value=DummyAlignmentFile(),
        ),
        patch("codfreq.sam2codfreq.JsonProgress", return_value=json_progress),
        patch("codfreq.sam2codfreq.ProcessPoolExecutor", DummyExecutor),
        patch(
            "codfreq.sam2codfreq.chunked_samfile", return_value=[(0, 0)]
        ),
        patch(
            "codfreq.sam2codfreq.codonalign_consensus",
            side_effect=lambda *a: a[:2],
        ),
    ):
        s2c.sam2codfreq(
            "file.bam",
            ref,
            cast(list, fragments),
            workers=1,
            log_format="json",
        )

    json_progress.update.assert_called_once_with(0)
    json_progress.close.assert_called_once()
