import unittest
from codfreq import posnas


class TestPosnas(unittest.TestCase):

    def test_merge_posnas(self) -> None:
        posnas1 = posnas.PosNA(10, 0, ord('A'))
        posnas2 = posnas.PosNA(10, 0, ord('C'))
        expected = posnas.PosNA(10, 0, ord('M'))
        self.assertEqual(posnas.merge_posnas(posnas1, posnas2), expected)

        posnas1 = posnas.PosNA(10, 0, ord('A'))
        posnas2 = posnas.PosNA(10, 0, ord('S'))
        expected = posnas.PosNA(10, 0, ord('V'))
        self.assertEqual(posnas.merge_posnas(posnas1, posnas2), expected)

        posnas1 = posnas.PosNA(10, 0, ord('A'))
        posnas2 = posnas.PosNA(10, 0, ord('A'))
        expected = posnas.PosNA(10, 0, ord('A'))
        self.assertEqual(posnas.merge_posnas(posnas1, posnas2), expected)

        posnas1 = posnas.PosNA(10, 0, ord('-'))
        posnas2 = posnas.PosNA(10, 0, ord('A'))
        expected = posnas.PosNA(10, 0, ord('-'))
        self.assertEqual(posnas.merge_posnas(posnas1, posnas2), expected)

        posnas1 = posnas.PosNA(10, 0, ord('W'))
        posnas2 = posnas.PosNA(10, 0, ord('M'))
        expected = posnas.PosNA(10, 0, ord('H'))
        self.assertEqual(posnas.merge_posnas(posnas1, posnas2), expected)
