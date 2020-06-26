import sys
import unittest

sys.path.append(".")
from reademptionlib.polyaclipper import PolyAClipper


class TestPolyAClipper(unittest.TestCase):
    def setUp(self):
        self.poly_a_clipper = PolyAClipper()

    def test_clip_poly_a_strech_no_change(self):
        """If there is no poly a strech the sequence is not changed."""
        test_seq = "ATAGTAGGAGATTTAGACCAGATGACGATGACACACAATGGTTAGGTACAGATAG"
        result_seq = test_seq
        self.assertEqual(
            self.poly_a_clipper.clip_poly_a_strech(test_seq), result_seq
        )

    def test_clip_poly_a_strech_empt(self):
        """If there the sequence is empty the sequence is not changed."""
        test_seq = ""
        result_seq = test_seq
        self.assertEqual(
            self.poly_a_clipper.clip_poly_a_strech(test_seq), result_seq
        )

    def test_clip_poly_a_strech_terminal_10_a(self):
        """Clipp terminal A strech if it is 10 A."""
        test_seq = "ATAGTAGGAGATTTAGACCAGATGACGATGACACAAAAAAAAAAA"
        result_seq = "ATAGTAGGAGATTTAGACCAGATGACGATGACAC"
        self.assertEqual(
            self.poly_a_clipper.clip_poly_a_strech(test_seq), result_seq
        )

    def test_clip_poly_a_strech_internal_10_a(self):
        """Clip before a 10 A long internal strech."""
        test_seq = "ATAGTAGGAGATTTAGACCAGATGACGATGACACAAAAAAAAAATTTAGACGACG"
        result_seq = "ATAGTAGGAGATTTAGACCAGATGACGATGACAC"
        self.assertEqual(
            self.poly_a_clipper.clip_poly_a_strech(test_seq), result_seq
        )

    def test_clip_poly_a_strech_terminal_09_a(self):
        """If there less than 10 A don't clip."""
        test_seq = "ATAGTAGGAGATTTAGACCAGATGACGATGACACAAAAAAAAAA"
        result_seq = test_seq
        self.assertEqual(
            self.poly_a_clipper.clip_poly_a_strech(test_seq), result_seq
        )

    def test_clip_poly_a_strecht_internal_(self):
        """If there less than 10 A don't clip."""
        test_seq = "ATAGTAGGAGATTTAGACCAGATGACGATGACACAAAAAAAAATTTAGACGACG"
        result_seq = test_seq
        self.assertEqual(
            self.poly_a_clipper.clip_poly_a_strech(test_seq), result_seq
        )

    def test_clip_poly_a_strecht_internal_09_a(self):
        """Test all A string"""
        test_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        result_seq = ""
        self.assertEqual(
            self.poly_a_clipper.clip_poly_a_strech(test_seq), result_seq
        )

    def test_aaaa_starting_substrings(self):
        test_seq = "TTTAAAATTTTTTTTAAAACCCCCCCCCCAAAAC"
        self.assertEqual(
            list(self.poly_a_clipper._aaaa_starting_substrings(test_seq, 11)),
            [["AAAATTTTTTT", 3], ["AAAACCCCCCC", 15]],
        )

    def test_remove_3_prime_a_no_change(self):
        """If there are no terminal As, there is no clipping."""
        test_seq = "AAAAATTTTCCGCCCGGGAAATTTT"
        result_seq = test_seq
        self.assertEqual(
            self.poly_a_clipper.remove_3_prime_a(test_seq), result_seq
        )

    def test_remove_3_prime_a_one_a(self):
        """Remove terminal A"""
        test_seq = "AAAAATTTTCCGCCCGGGAAATTTTA"
        result_seq = "AAAAATTTTCCGCCCGGGAAATTTT"
        self.assertEqual(
            self.poly_a_clipper.remove_3_prime_a(test_seq), result_seq
        )

    def test_remove_3_prime_a_multiple_as(self):
        """Remove terminal stretch of multiple As"""
        test_seq = "AAAAATTTTCCGCCCGGGAAATTTTAAAAAA"
        result_seq = "AAAAATTTTCCGCCCGGGAAATTTT"
        self.assertEqual(
            self.poly_a_clipper.remove_3_prime_a(test_seq), result_seq
        )


if __name__ == "__main__":
    unittest.main()
