import unittest
import sys
sys.path.append('.')
sys.path.append('..')
from poly_a_clipper import PolyAClipper

class OptionsMock:
    tmp = 2
    output_file = None

class TestPolyAClipper(unittest.TestCase):
    """ 
    """

    def setUp(self):
        args = ["humpy-dumpy.fa"]
        options = OptionsMock()
        self.poly_a_clipper = PolyAClipper(args, options)

    def test__test_and_clip_1(self):
        """ Test _test_and_clip

        No reason to clip.
        """
        input_seq = "TTTAAATTTCACATATAGATAGTATAGGA"
        expected_output_seq = input_seq
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))

    def test__test_and_clip_2(self):
        """Test _test_and_clip

        AAAA but nothing else.
        """
        #               1234
        #               ====
        input_seq = "TTTAAAACACATATGTATAGGA"
        expected_output_seq = input_seq
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))

    def test__test_and_clip_3(self):
        """Test _test_and_clip

        AAAAAAAA but nothing else.
        """
        #               12345678
        #               --------
        input_seq = "TTTAAAAAAAACATATAAATAGATAGTATAGGA"
        expected_output_seq = input_seq
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))

    def test__test_and_clip_4(self):
        """Test _test_and_clip

        AAAA and AAAAAAAAA without any mismatch, insertion or
        deletion.
        """
        #                         123412345678
        #                         ====--------
        input_seq = "ATCTATGGCGCGCAAAAAAAAAAAATGATGAAAAGGTC"
        expected_output_seq = "ATCTATGGCGCGC"
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))

    def test__test_and_clip_5(self):
        """Test _test_and_clip
        AAAA and AAAAAAAAA (A x 8) with 1 mismatch in the middle
        """
        #                         123412345678
        #                         ====---X----
        input_seq = "ATCTATGGCGCGCAAAAAAATAAAATGATGAACAAGGTC"
        expected_output_seq = "ATCTATGGCGCGC"
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))

    def test__test_and_clip_6(self):
        """Test _test_and_clip

        AAAA and AAAAAAAAA (A x 8) with 1 mismatch in the end
        """
        #                         123412345678
        #                         ====-------X
        input_seq = "ATCTATGGCGCGCAAAAAAAAAAACTGATGAAAAGGTC"
        expected_output_seq = "ATCTATGGCGCGC"
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))

    def test__test_and_clip_7(self):
        """Test _test_and_clip

        AAAA and AAAAAAAAA (8 characters) with 2 mismatches.
        """
        #                         123412345678
        #                         ====---X--X-
        input_seq = "ATCTATGGCGCGCAAAAAAATAACATGATGAACAAGGTC"
        expected_output_seq = input_seq
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))

    def test__test_and_clip_8(self):
        """Test _test_and_clip

        AAAA and AAAAAAAAA (9 characters) with 1 mismatch.
        """
        #                             1234123456789
        #                             ====-X-------
        input_seq = "GCGGGATATAGTACGACAAAAAGAAAAAAAGCCCTTGCTAT "
        expected_output_seq = "GCGGGATATAGTACGAC"
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))

    def test__test_and_clip_9(self):
        """Test _test_and_clip

        AAAA and AAAAAAAAAA (9 characters) with 2 mismatches.
        """
        #                               1234123456789
        #                               ====--X--X---
        input_seq = "GCGGGATATAGTACGACAAAAAATAAGAAAGCCCTTGCTAT"
        expected_output_seq = input_seq
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))
        

    def test__test_and_clip_10(self):
        """Test _test_and_clip

        AAAA and AAAAAAA (7 characters) with 1 mismatch.
        """
        #                          12341234567
        #                          ====----X-- 
        input_seq = "GAGTCATACCGTGCAAAAAAAACAAGCGTTTGAGATT"
        expected_output_seq = "GAGTCATACCGTGC"
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))


    def test__test_and_clip_11(self):
        """Test _test_and_clip

        AAAA and AAAAAAA (7 characters) with 2 mismatch.
        """
        #                          12341234567
        #                          ====--X-X-- 
        input_seq = "GAGTCATACCGTGCAAAAAAGACAAGCGTTTGAGATT"
        expected_output_seq = input_seq
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))



    def test__test_and_clip_12(self):
        """Test _test_and_clip

        AAAA and AAAAAAAAA (A x 8) with one AAAA in before this pattern.
        """
        #                         1234     123412345678
        #                         ====     ====--------
        input_seq = "AGAGCAGCGCATTAAAAGCGCGAAAAAAAAAAAACGCGTGA"
        expected_output_seq = "AGAGCAGCGCATTAAAAGCGCG"
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))

    def test__test_and_clip_13(self):
        """Test _test_and_clip

        AAAA and AAAAAAAAA (A x 8) with one AAAA in behind this pattern.
        """
        #                              123412345678      1234
        #                              ====--------      ====
        input_seq = "AGAGCAGCGCATTGCGCGAAAAAAAAAAAACGCGTGAAAATGCT"
        expected_output_seq = "AGAGCAGCGCATTGCGCG"
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))

    def test__test_and_clip_14(self):
        """Test _test_and_clip

        AAAA and AAAAAAAAA (A x 8) at the very beginning of the
        sequence.
        """
        #            123412345678
        #            ====--------
        input_seq = "AAAAAAAAAAAAGCGGGGGGCGGG"
        expected_output_seq = "A"
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))

    def test__test_and_clip_15(self):
        """Test _test_and_clip

        AAAA very close to the end of the string.
        """
        #                  123412345
        #                  ====-----
        input_seq = "TTTGTGAAAAAAGGG"
        expected_output_seq = "TTTGTGAAAAAAGGG"
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))


    def test__test_and_clip_16(self):
        """Test _test_and_clip

        AAAA and AAAAAAAAA (A x 8) at the very beginning of the
        sequence.
        """
        #                          123412345678912345
        #                          ====--------------
        input_seq = "AATGTGTGTGTGTGAAAAAAAAAAAAAAAAAAGCGGGGGGCGGG"
        expected_output_seq = "AATGTGTGTGTGTG"
        self.assertEqual(
            expected_output_seq, 
            self.poly_a_clipper._test_and_clip(input_seq))

    def test__aaaa_starting_substrings_1(self):
        """Test _aaaa_starting_substrings
        
        The sequence has the AAAA at the very beginnig.

        """
        seq = "AAAATGTGTTTGTGGTTTTGGTG"
        length = 11
        seqs = [seq for seq in self.poly_a_clipper._aaaa_starting_substrings(
                seq, length)]
        self.assertEqual([["AAAATGTGTTT", 0]], seqs)

    def test__aaaa_starting_substrings_2(self):
        """Test _aaaa_starting_substrings

        The sequence has the AAAA at the very beginnig and directly
        another A following. This should lead to to two subsequences.
        """
        seq = "AAAAATGTGTTTGTGGTTTTGGTG"
        length = 11
        seqs = [seq for seq in self.poly_a_clipper._aaaa_starting_substrings(
                seq, length)]
        self.assertEqual([["AAAAATGTGTT", 0], 
                          ["AAAATGTGTTT", 1]], seqs)

    def test__aaaa_starting_substrings_3(self):
        """Test _aaaa_starting_substrings

        The sequence has the AAAA at the very beginnig and directly
        another A following. This should lead to to two
        subsequences. Additionally a AAA occurs in the middel of the
        sequence leading to a third subsequence.
        """
        seq = "AAAAATGTGTTTGTGGTTTTGGTGAAAAGGGTGAGATTGATGGTAAG"
        length = 11
        seqs = [seq for seq in self.poly_a_clipper._aaaa_starting_substrings(
                seq, length)]
        self.assertEqual([["AAAAATGTGTT", 0], ["AAAATGTGTTT", 1], 
                          ["AAAAGGGTGAG", 24]], seqs)

    def test__aaaa_starting_substrings_4(self):
        """Test _aaaa_starting_substrings
        
        No AAAA in the sequence => not hit.
        """
        seq = "TGTGTTTGTGGTTTTGGTGGTGAGATTGATGGTAAG"
        length = 11
        seqs = [seq for seq in self.poly_a_clipper._aaaa_starting_substrings(
                seq, length)]
        self.assertEqual([], seqs)

    def test__aaaa_starting_substrings_5(self):
        """Test _aaaa_starting_substrings
        
        The AAA occure too close at the end of the sequence => no hit.
        """
        seq = "TGTGTTTGTGGTTTTGGTGGTGAGATTGATGGTAAGAAAATTG"
        length = 11
        seqs = [seq for seq in self.poly_a_clipper._aaaa_starting_substrings(
                seq, length)]
        self.assertEqual([], seqs)

if __name__ == "__main__":
    unittest.main()
