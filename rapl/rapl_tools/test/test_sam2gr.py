import unittest
import sys
sys.path.append('.')
sys.path.append('..')
from sam2gr import Sam2Gr

class ArgsMock_1(object):
    sam_file = "humpy-dumpy.sam"
    normalize = False
    output_prefix = None
    normalization_factor = None
    normalization_multiplier = None
    normalize_by_reads = False
    normalize_by_nucleotides = False
    starts_only = False
    mapping_target = None

class TestSam2Gr_1(unittest.TestCase):
    """Test Sam2Gr"""

    def setUp(self):
        args = ArgsMock_1()
        self.sam2gr = Sam2Gr(
            args.sam_file, args.mapping_target, args.normalize_by_reads,
            args.normalize_by_nucleotides, args.normalization_factor,
            args.normalization_multiplier, args.output_prefix, args.starts_only)

    def test__extension_for_intensities_list_1(self):
        """Test: The current intensity list is 4 elements long. 

        The last position of a new entrie is 7. So we need to extend
        the list by an extension list that is 3 elements long.

        """
        intensities = [1, 2, 3 , 2]
        ext_intensities = self.sam2gr._extension_for_intensities_list(
            intensities, 7)
        self.assertEqual([0, 0, 0], ext_intensities)

    def test__extension_for_intensities_list_2(self):
        """Test: The current intensity list is 4 elements long. 

        The last position of a new entrie is 5. So we need to extend
        the list by an extension list that is 1 element long.

        """
        intensities = [1, 2, 3 , 2]
        ext_intensities = self.sam2gr._extension_for_intensities_list(
            intensities, 5)
        self.assertEqual([0], ext_intensities)
        
    def test__add_value_to_coordinates_1(self):
        """Test the addition of 1 to some positions."""
        intensities =  [1, 1, 1, 1, 1]
        start = 2
        end = 4
        strand = "+"
        number_of_hits = 1
        self.sam2gr._add_value_to_coordinates(
            intensities,  start, end, strand, number_of_hits)
        self.assertEqual([1, 2, 2, 2, 1], intensities)

    def test__add_value_to_coordinates_2(self):
        """Test the addition of 1 to some element and the extension of
        the list.

        """
        intensities =  [1, 1, 1, 1, 1]
        start = 3
        end = 5
        strand = "+"
        number_of_hits = 1
        self.sam2gr._add_value_to_coordinates(
            intensities, start, end, strand, number_of_hits)
        self.assertEqual([1, 1, 2, 2, 2], intensities)

    def test__add_value_to_coordinates_3(self):
        """Test the addition of 1 to some element and the extension of
        the list.

        """
        intensities =  [1, 1, 1, 1, 1]
        start = 6
        end = 8
        strand = "+"
        number_of_hits = 1
        self.sam2gr._add_value_to_coordinates(
            intensities, start, end, strand, number_of_hits)
        self.assertEqual([1, 1, 1, 1, 1, 1, 1, 1], intensities)

    def test__add_value_to_coordinates_2(self):
        """Test the substract of 1 to some positions (minus strand)."""
        intensities =  [-1, -1, -1, -1, -1]
        start = 2
        end = 4
        strand = "-"
        number_of_hits = 1
        self.sam2gr._add_value_to_coordinates(
            intensities, start, end, strand, number_of_hits)
        self.assertEqual([-1, -2.0, -2.0, -2.0, -1], intensities)

    def test__normalized_intensity_1(self):
        """Test the normalization of the intensities by the default value."""
        intensity = 20
        self.sam2gr.final_normalization_factor = 1
        self.assertEqual(20, self.sam2gr._normalized_intensity(intensity))

    def test__normalized_intensity_2(self):
        """Test the normalization of the intensities by a given value."""
        intensity = 20
        self.sam2gr.final_normalization_factor = 0.5
        self.assertEqual(10, self.sam2gr._normalized_intensity(intensity))

    def test__generate_suffix_1(self):
        """Test the output file suffix for runs without normalization"""
        self.assertEqual("_plus.gr", self.sam2gr._generate_suffix("plus"))
        self.assertEqual("_minus.gr", 
                         self.sam2gr._generate_suffix("minus"))

    def test__process_entry_plus_strand_adds(self):
        self.sam2gr.intensities_plus =  []
        self.sam2gr.intensities_minus =  []
        sam_entry = {"query" : "IPAR1_0028_FC:5:106:17555:7748#0", "flag" : 0,
                     "start" : 3, "number_of_hits": "NH:i:1", "cigar" : "7M"}
        self.sam2gr._process_entry(sam_entry)
        self.assertEqual([0, 0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 
                         self.sam2gr.intensities_plus)
        self.assertEqual([], self.sam2gr.intensities_minus)

    def test__process_entry_minus_start_adds(self):
        self.sam2gr.intensities_plus =  []
        self.sam2gr.intensities_minus =  []
        sam_entry = {"query" : "IPAR1_0028_FC:5:106:17555:7748#0", "flag" : 16, 
                     "start" : 2, "number_of_hits": "NH:i:1", "cigar" : "5M"}
        self.sam2gr._process_entry(sam_entry)
        self.assertEqual([], self.sam2gr.intensities_plus)
        self.assertEqual([0, -1.0, -1.0, -1.0, -1.0, -1.0], 
                         self.sam2gr.intensities_minus)

    def test__process_entry_multiple_mapped_read(self):
        self.sam2gr.intensities_plus =  []
        self.sam2gr.intensities_minus =  []
        sam_entry = {"query" : "IPAR1_0028_FC:5:106:17555:7748#0", "flag" : 0,
                     "start" : 3, "number_of_hits": "NH:i:5", "cigar" : "7M"}
        self.sam2gr._process_entry(sam_entry)
        self.assertEqual([0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2], 
                         self.sam2gr.intensities_plus)
        self.assertEqual([], self.sam2gr.intensities_minus)        
        
class ArgMock_2(object):
    sam_file = "humpy-dumpy.sam"
    normalize = True
    normalization_factor = 50
    normalization_multiplier = 1000
    output_prefix = None
    normalize_by_reads = False
    normalize_by_nucleotides = False
    starts_only = False
    mapping_target = None

class TestSam2Gr_2(unittest.TestCase):
    """Test Sam2Gr."""

    def setUp(self):
        args = ArgMock_2()
        self.sam2gr = Sam2Gr(
            args.sam_file, args.mapping_target, args.normalize_by_reads,
            args.normalize_by_nucleotides, args.normalization_factor,
            args.normalization_multiplier, args.output_prefix, 
            args.starts_only)

    def test__set_final_normalization_factor_3(self):
        """Test the calculation of the final normalization factor for
        the case that the normalization should be done based on the
        given normalization factor and multiplier.
        """
        self.sam2gr.normalize = True
        self.sam2gr._set_final_normalization_factor()
        self.assertEqual(50, self.sam2gr.normalization_factor)
        self.assertEqual(1000, self.sam2gr.normalization_multiplier)
        self.assertEqual(20.0, self.sam2gr.final_normalization_factor)

class ArgMock_3(object):
    sam_file = "humpy-dumpy.sam"
    normalize = True
    normalization_factor = None
    normalization_multiplier = None
    output_prefix = None
    normalize_by_reads = False
    normalize_by_nucleotides = False
    starts_only = False
    mapping_target = None
    normalization_multiplier = 1000.0

class TestSam2Gr_3(unittest.TestCase):
    """Test Sam2Gr with normalization by read and default
    multiplier.

    """
    
    def setUp(self):
        args = ArgMock_3()
        self.sam2gr = Sam2Gr(
            args.sam_file, args.mapping_target, args.normalize_by_reads,
            args.normalize_by_nucleotides, args.normalization_factor,
            args.normalization_multiplier, args.output_prefix, 
            args.starts_only)
        
    def test__set_final_normalization_factor_1(self):
        """Test the calculation of the final normalization factor for
        the case that the normalization should be done based on the
        number of reads.
        """
        self.sam2gr.normalize = True
        self.sam2gr.reads_and_countings = {
            "x" : 1, "y" : 1, "z" : 1, "p" : 1}
        self.sam2gr.normalization_factor = 4.0
        self.sam2gr._set_final_normalization_factor()
        self.assertEqual(1000.0, self.sam2gr.normalization_multiplier)
        self.assertEqual(4, len(self.sam2gr.reads_and_countings))
        self.assertEqual(250.0, self.sam2gr.final_normalization_factor)

    def test__generate_suffix_1(self):
        """Test the output file suffix for runs with normalization and
        calculated normalization factor.
        """
        self.sam2gr.normalization_factor = 3.0
        self.sam2gr.normalization_factor_for_suffix = 3.0
        self.sam2gr.normalization_multiplier_for_suffix = 1000.0
        self.assertEqual("_normalized_by_3.0_multiplied_by_1000.0_plus.gr", 
                         self.sam2gr._generate_suffix("plus"))
        self.assertEqual("_normalized_by_3.0_multiplied_by_1000.0_minus.gr", 
                         self.sam2gr._generate_suffix("minus"))
        
if __name__ == "__main__":
    unittest.main()
