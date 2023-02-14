import os
from typing import Tuple, Optional

import unittest
from codfreq.posnas import join_posnas, PosNA
from codfreq.segfreq import SegFreq, remove_first_n_pos, remove_last_n_pos

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
SARS2_SEGFREQ = os.path.join(DATA_DIR, 'sars2-test.segfreq')
SARS2_FASTA = os.path.join(DATA_DIR, 'sars2-test.fasta')


class TestSegFreq(unittest.TestCase):

    def setUp(self) -> None:
        with open(SARS2_SEGFREQ, encoding='UTF-8-sig') as fh:
            self.segfreq = SegFreq.load(fh)

    def test_get_frequency(self) -> None:
        self.assertEqual(
            self.segfreq.get_frequency([21563]),
            {b'ATG': 1377, b'ACG': 2, b'GTG': 2}
        )

        self.assertEqual(
            self.segfreq.get_frequency([21563 + 569 * 3]),
            {b'GAT': 2833, b'GCT': 122, b'GGT': 13, b'GTT': 1,
             b'TAT': 1, b'AAT': 2, b'CAA': 1, b'GAA': 1, b'GAC': 7}
        )

    def test_get_pos_nas(self) -> None:
        self.assertEqual(
            self.segfreq.get_pos_nas(21563),
            {b'A': 1411, b'G': 2}
        )

    def test_get_consensus(self) -> None:
        self.assertEqual(self.segfreq.segment_size, 12)
        self.assertEqual(self.segfreq.segment_step, 4)
        with open(SARS2_FASTA) as fh:
            seq = fh.read().strip()
            self.assertEqual(
                join_posnas(self.segfreq.get_consensus(21563, 25381)),
                seq
            )

    def test_get_patterns(self) -> None:
        # codon
        self.assertEqual(
            {
                (posnas[0].pos, join_posnas(posnas)): count
                for posnas, count in
                self.segfreq.get_patterns(
                    21563, 21563 + 2, top_n_seeds=10).items()
            },
            {
                (21563, 'ATG'): (1377, .9686),
                (21563, 'AT'): (26, .0182),
                (21564, 'TG'): (6, .0042),
                (21563, 'A'): (6, .0042),
                (21563, 'GTG'): (2, .0014),
                (21563, 'ACG'): (2, .0014)
            }
        )
        # single segment
        self.assertEqual(
            {
                (posnas[0].pos, join_posnas(posnas)): count
                for posnas, count in
                self.segfreq.get_patterns(
                    21563, 21563 + 4, top_n_seeds=10).items()
            },
            {
                (21563, 'ATGTT'): (1354, .9525),
                (21563, 'AT'): (26, .0182),
                (21563, 'ATG'): (10, .007),
                (21563, 'ATGTC'): (7, .0049),
                (21564, 'TGTT'): (6, .0042),
                (21563, 'A'): (6, .0042),
                (21563, 'ATGCT'): (2, .0014),
                (21563, 'GTGTT'): (2, .0014),
                (21563, 'ATGT'): (2, .0014),
                (21563, 'ACGTT'): (2, .0014)
            }
        )

        # merge two segments
        self.assertEqual(
            {
                (posnas[0].pos, join_posnas(posnas)): count
                for posnas, count in
                self.segfreq.get_patterns(
                    21563, 21563 + 10, top_n_seeds=5).items()
            },
            {
                (21563, 'ATGTTTGTTTT'): (1320, .9287),
                (21563, 'AT'): (26, .0182),
                (21564, 'TGTTTGTTTT'): (6, .0042),
                (21563, 'GTGTTTGTTTT'): (2, .0014),
                (21563, 'ACGTTTGTTTT'): (1, .0007)
            }
        )

        # high variation segments
        self.assertEqual(
            {
                (posnas[0].pos, join_posnas(posnas)): count
                for posnas, count in
                self.segfreq.get_patterns(
                    21563 + 568 * 3,
                    21563 + 568 * 3 + 24,
                    top_n_seeds=5
                ).items()
            },
            {
                (23267, 'ATTGATGACACTACTGATGCTGTCC'): (2563, .7837),
                (23267, 'ATTG'): (173, .0511),
                (23267, 'ATTGCTGACACTACTGATGCTGTCC'): (114, .0342),
                (23279, 'ACTGATGCTGTCC'): (85, .0264),
                (23272, 'TGACACTACTGATGCTGTCC'): (13, .0192)
            }
        )

    def test_remove_first_n_pos(self) -> None:
        segments: Tuple[Optional[PosNA], ...] = (
            PosNA(7, 0, ord('A')),
            PosNA(7, 1, ord('A')),
            PosNA(8, 0, ord('A')),
            PosNA(8, 1, ord('A')),
            PosNA(9, 0, ord('A')),
            PosNA(9, 1, ord('A')),
            PosNA(10, 0, ord('A')),
            PosNA(10, 1, ord('A'))
        )
        self.assertEqual(
            remove_first_n_pos(segments, 2),
            segments[4:]
        )
        segments = (None,) * 2 + segments
        self.assertEqual(
            remove_first_n_pos(segments, 4),
            segments[6:]
        )

    def test_remove_last_n_pos(self) -> None:
        segments: Tuple[Optional[PosNA], ...] = (
            PosNA(7, 0, ord('A')),
            PosNA(7, 1, ord('A')),
            PosNA(8, 0, ord('A')),
            PosNA(8, 1, ord('A')),
            PosNA(9, 0, ord('A')),
            PosNA(9, 1, ord('A')),
            PosNA(10, 0, ord('A')),
            PosNA(10, 1, ord('A'))
        )
        self.assertEqual(
            remove_last_n_pos(segments, 2),
            segments[:4]
        )
        segments = segments + (None,) * 2
        self.assertEqual(
            remove_last_n_pos(segments, 4),
            segments[:4]
        )
