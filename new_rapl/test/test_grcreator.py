import sys
from io import StringIO
import unittest
sys.path.append(".")
from libs.grcreator import GRCreator, GRFileBuilder

class TestGRFileBuilder(unittest.TestCase):

    def setUp(self):
        self.gr_file_builder = GRFileBuilder(
            "a_sam_file", "a_ref_genome", "output_plus", "output_minus",
            normalization_value=1, multiplier=2)
        
    def test_calc_raw_coverages(self):
        self.gr_file_builder._sam_entries = mock_sam_entries
        self.assertTupleEqual(
            ([1.0, 1.0, 2.0, 2.0, 2.0, 1.0, 1.0], 
             [-1.0, -2.0, -2.0, -2.0, -2.0, -1.0]), 
            self.gr_file_builder._calc_raw_coverages("a_handle_spaceholder"))

    def test_add_coverage(self):
        # Add to empty coverage list
        sam_entry = MockSamEntry(1, 2, "a_ref_genome", "+", 1)
        coverages = []
        self.gr_file_builder._add_coverage(sam_entry, coverages)
        self.assertListEqual([1.0, 1.0], coverages)
        # Add to filled coverage list - plus strand
        sam_entry = MockSamEntry(1, 2, "a_ref_genome", "+", 1)
        coverages = [1.0, 1.0]
        self.gr_file_builder._add_coverage(sam_entry, coverages)
        self.assertListEqual([2.0, 2.0], coverages)
        # Add to filled coverage list - minus strand
        sam_entry = MockSamEntry(1, 2, "a_ref_genome", "-", 1)
        coverages = [1.0, 1.0]
        self.gr_file_builder._add_coverage(sam_entry, coverages)
        self.assertListEqual([2.0, 2.0], coverages)
        # Add to filled coverage list with more multi hit entry
        sam_entry = MockSamEntry(1, 2, "a_ref_genome", "+", 4)
        coverages = [1.0, 1.0]
        self.gr_file_builder._add_coverage(sam_entry, coverages)
        self.assertListEqual([1.25, 1.25], coverages)
        # Add inculing list extension
        sam_entry = MockSamEntry(1, 4, "a_ref_genome", "+", 1)
        coverages = [1.0, 1.0]
        self.gr_file_builder._add_coverage(sam_entry, coverages)
        self.assertListEqual([2.0, 2.0, 1.0, 1.0], coverages)
        # Add values in the middle of the list
        sam_entry = MockSamEntry(2, 3, "a_ref_genome", "+", 1)
        coverages = [1.0, 1.0, 1.0, 1.0]
        self.gr_file_builder._add_coverage(sam_entry, coverages)
        self.assertListEqual([1.0, 2.0, 2.0, 1.0], coverages)

    def test__normalize_and_multiply(self):
        # Test multiplier
        self.gr_file_builder.multiplier = 2
        self.gr_file_builder.normalization_value = 1
        coverages = [1.0, 2.0, 4.0, 8.0]
        self.assertListEqual(
            [2.0, 4.0, 8.0, 16.0], 
            self.gr_file_builder._normalize_and_multiply(coverages))
        # Test normalization value
        self.gr_file_builder.multiplier = 1
        self.gr_file_builder.normalization_value = 2
        coverages = [1.0, 2.0, 4.0, 8.0]
        self.assertListEqual(
            [0.5, 1.0, 2.0, 4.0], 
            self.gr_file_builder._normalize_and_multiply(coverages))
        # Test multiplier and  normalization value
        self.gr_file_builder.multiplier = 10
        self.gr_file_builder.normalization_value = 2
        coverages = [1.0, 2.0, 4.0, 8.0]
        self.assertListEqual(
            [5.0, 10.0, 20.0, 40.0], 
            self.gr_file_builder._normalize_and_multiply(coverages))
        
    def test_extend_coverages(self):
        coverages = [1.0, 1.0]
        self.gr_file_builder._extend_coverages(coverages, 4)
        self.assertListEqual([1.0, 1.0, 0.0, 0.0], coverages)

    def test_norm_or_multi_needed(self):
        # Should be False per default
        gr_file_builder = GRFileBuilder(
            "a_sam_file", "a_ref_genome", "output_plus", "output_minus")
        self.assertEqual(1, gr_file_builder.multiplier)
        self.assertEqual(1, gr_file_builder.normalization_value)
        self.assertFalse(gr_file_builder._norm_or_multi_needed())
        # If normalization_value != 1 it should be needed
        gr_file_builder.multiplier = 1
        gr_file_builder.normalization_value = 2.0
        self.assertTrue(gr_file_builder._norm_or_multi_needed())
        # If multiplier != 1 it should be needed
        gr_file_builder.multiplier = 2.0
        gr_file_builder.normalization_value = 1
        self.assertTrue(gr_file_builder._norm_or_multi_needed())
        # If both value are != 0 it should be needed
        gr_file_builder.multiplier = 2.0
        gr_file_builder.normalization_value = 1
        self.assertTrue(gr_file_builder._norm_or_multi_needed())        

class MockSamEntry(object):
    
    def __init__(self, start, end, reference, strand, number_of_hits_as_int):
        self.start = start
        self.end = end
        self.reference = reference
        self.strand = strand
        self.number_of_hits_as_int = number_of_hits_as_int

def mock_sam_entries(sam_fh):
    for entry in [
        MockSamEntry(1, 5, "a_ref_genome", "+", 1),
        MockSamEntry(1, 10, "another_ref_genome", "+", 1),
        MockSamEntry(3, 7, "a_ref_genome", "+", 1),
        MockSamEntry(1, 5, "a_ref_genome", "-", 1),
        MockSamEntry(2, 6, "a_ref_genome", "-", 1)]:
        yield(entry)

if __name__ == "__main__":
    unittest.main()
