from collections import Counter
from typing import cast

from unittest.mock import MagicMock, patch

import codfreq.sam2codfreq as s2c
from codfreq.codfreq_types import (
    DerivedFragmentConfig,
    MainFragmentConfig,
    Profile,
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
    profile = {
        "fragmentConfig": [
            {"fragmentName": "refA", "refSequence": "AAA"},
            {
                "fragmentName": "fragA",
                "fromFragment": "refA",
                "geneName": "geneX",
                "refRanges": [(1, 3)],
                "codonAlignment": [{
                    "relRefStart": 1,
                    "relRefEnd": 3,
                    "windowSize": 3,
                    "minGapDistance": 5,
                    "relGapPlacementScore": "0:0-1=1",
                }],
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
    refs, lookup = s2c.get_ref_fragments(profile)
    assert refs[0][2][0]["codonAlignment"][0]["relRefEnd"] == 3
    assert refs[0][2][0]["codonAlignment"][0]["minGapDistance"] == 5
    assert (
        refs[0][2][0]["codonAlignment"][0]["relGapPlacementScore"] == "0:0-1=1"
    )
    assert refs[0][2][1]["codonAlignment"] is False
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


def test_sam2codfreq() -> None:
    """Combine chunk results into codon counters."""

    alignment_mock = MagicMock(mapped=3)
    alignment_mock.__enter__.return_value = alignment_mock
    alignment_mock.__exit__.return_value = False

    pbar_mock = MagicMock()

    executor_inst = MagicMock()
    executor_inst.__enter__.return_value = executor_inst
    executor_inst.__exit__.return_value = False
    executor_inst.map.side_effect = (
        lambda fn, *iterables: [fn(*args) for args in zip(*iterables)]
    )

    with patch(
        "codfreq.sam2codfreq.pysam.AlignmentFile",
        return_value=alignment_mock,
    ), patch("codfreq.sam2codfreq.tqdm", return_value=pbar_mock), patch(
        "codfreq.sam2codfreq.chunked_samfile", return_value=[(0, 1)]
    ), patch(
        "codfreq.sam2codfreq.ProcessPoolExecutor", return_value=executor_inst
    ), patch.object(
        s2c,
        "sam2codfreq_between",
        return_value=(
            Counter({("fragA", 1, "AAA"): 2}),
            Counter({("fragA", 1, "AAA"): 40}),
            2,
        ),
    ), patch.object(
        s2c,
        "codonalign_consensus",
        side_effect=lambda stat, qual, ref, frags: (stat, qual),
    ):
        ref = cast(
            MainFragmentConfig, {"fragmentName": "refA", "refSequence": "AAA"}
        )
        fragments = [
            cast(
                DerivedFragmentConfig,
                {
                    "fragmentName": "fragA",
                    "fromFragment": "refA",
                    "refRanges": [(1, 3)],
                },
            )
        ]
        stat_by_fragpos, qual_by_fragpos = s2c.sam2codfreq(
            "file.bam", ref, fragments, workers=1
        )

    assert stat_by_fragpos == {("fragA", 1): Counter({"AAA": 2})}
    assert qual_by_fragpos == {("fragA", 1): Counter({"AAA": 40})}


def test_sam2codfreq_json_logging() -> None:
    """JSON log format uses JsonProgress."""

    alignment_mock = MagicMock(mapped=2)
    alignment_mock.__enter__.return_value = alignment_mock
    alignment_mock.__exit__.return_value = False

    pbar_mock = MagicMock()

    executor_inst = MagicMock()
    executor_inst.__enter__.return_value = executor_inst
    executor_inst.__exit__.return_value = False
    executor_inst.map.side_effect = (
        lambda fn, *iterables: [fn(*args) for args in zip(*iterables)]
    )

    with (
        patch(
            "codfreq.sam2codfreq.pysam.AlignmentFile",
            return_value=alignment_mock,
        ),
        patch("codfreq.sam2codfreq.JsonProgress", return_value=pbar_mock),
        patch("codfreq.sam2codfreq.chunked_samfile", return_value=[(0, 1)]),
        patch(
            "codfreq.sam2codfreq.ProcessPoolExecutor",
            return_value=executor_inst,
        ),
        patch.object(
            s2c,
            "sam2codfreq_between",
            return_value=(
                Counter({
                    ("fragA", 1, "AAA"): 1,
                }),
                Counter({
                    ("fragA", 1, "AAA"): 20,
                }),
                1,
            ),
        ),
        patch.object(
            s2c,
            "codonalign_consensus",
            side_effect=lambda stat, qual, ref, frags: (stat, qual),
        ),
    ):
        ref = cast(
            MainFragmentConfig,
            {"fragmentName": "refA", "refSequence": "AAA"},
        )
        fragments = [
            cast(
                DerivedFragmentConfig,
                {
                    "fragmentName": "fragA",
                    "fromFragment": "refA",
                    "refRanges": [(1, 3)],
                },
            )
        ]
        s2c.sam2codfreq(
            "file.bam", ref, fragments, workers=1, log_format="json"
        )

    pbar_mock.update.assert_called_once_with(1)
    pbar_mock.close.assert_called_once()


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
        profile = cast(
            Profile,
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
            },
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
