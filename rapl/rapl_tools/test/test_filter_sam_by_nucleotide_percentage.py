from io import StringIO
import sys
import unittest

sys.path.append('..')
sys.path.append('.')
from filter_sam_by_nucleotide_percentage import SamNuclPercFilter

class OptionsMock(object):
    output_file_prefix = None


class TestSamNucleotidePercFilter(unittest.TestCase):
    """
    """

    def setUp(self):
        args = ["boing", "A", 40]
        options = OptionsMock()
        self.sam_nucl_perc_filter = SamNuclPercFilter(
            args, options)

    def test__nucl_percentage_too_high_1(self):
        self.sam_nucl_perc_filter.max_percentage = 50
        self.sam_nucl_perc_filter.nucleotide = "A"
        sequence = "AAAAAAAAAAA"
        strand = "+"
        self.assertEqual(
            True, 
            self.sam_nucl_perc_filter._nucl_percentage_too_high(
                sequence, strand))

    def test__nucl_percentage_too_high_2(self):
        self.sam_nucl_perc_filter.max_percentage = 50
        self.sam_nucl_perc_filter.nucleotide = "A"
        sequence = "TTTTTTTTTT"
        strand = "+"
        self.assertEqual(
            False, 
            self.sam_nucl_perc_filter._nucl_percentage_too_high(
                sequence, strand))

    def test__nucl_percentage_too_high_2(self):
        self.sam_nucl_perc_filter.max_percentage = 50
        self.sam_nucl_perc_filter.nucleotide = "A"
        sequence = "TTTTTTTTTT"
        strand = "-"
        self.assertEqual(
            True, 
            self.sam_nucl_perc_filter._nucl_percentage_too_high(
                sequence, strand))

    def test__effective_nucleotide_1(self):
        self.sam_nucl_perc_filter.nucleotide = "A"
        strand = "+"
        self.assertEqual(
            "A", 
            self.sam_nucl_perc_filter._effective_nucleotide(strand))

    def test__effective_nucleotide_2(self):
        self.sam_nucl_perc_filter.nucleotide = "A"
        strand = "-"
        self.assertEqual(
            "T", 
            self.sam_nucl_perc_filter._effective_nucleotide(strand))

    def test__effective_nucleotide_3(self):
        self.sam_nucl_perc_filter.nucleotide = "C"
        strand = "+"
        self.assertEqual(
            "C", 
            self.sam_nucl_perc_filter._effective_nucleotide(strand))

    def test__effective_nucleotide_4(self):
        self.sam_nucl_perc_filter.nucleotide = "T"
        strand = "-"
        self.assertEqual(
            "A", 
            self.sam_nucl_perc_filter._effective_nucleotide(strand))

    def test__perc_in_sequence_1(self):
        nucleotide = "A"
        sequence = "AAAAAAA"
        self.assertEqual(
            100.0, 
            self.sam_nucl_perc_filter._perc_in_sequence(
                nucleotide, sequence))

    def test__perc_in_sequence_2(self):
        nucleotide = "A"
        sequence = "CTTGTGTCGT"
        self.assertEqual(
            0.0, 
            self.sam_nucl_perc_filter._perc_in_sequence(
                nucleotide, sequence))

    def test__perc_in_sequence_3(self):
        nucleotide = "C"
        sequence = "CCCCCTTTTT"
        self.assertEqual(
            50.0, 
            self.sam_nucl_perc_filter._perc_in_sequence(
                nucleotide, sequence))

    def test_process_sam_entry_with_passing_seq(self):
        self._prepare_for_processing_test()
        self.entry['sequence'] = "TGTGGGACATTGTGGATG"
        self.sam_nucl_perc_filter._process_sam_entry(self.entry)
        self.assertEqual(self.sam_nucl_perc_filter.counter_sufficient, 1)
        self.assertEqual(self.sam_nucl_perc_filter.counter_insufficient, 0)
        self.assertTrue(
            self.sam_nucl_perc_filter.output_fh_sufficient.getvalue().startswith(
                "test_seq"))
        self.assertEqual(
            self.sam_nucl_perc_filter.output_fh_insufficient.getvalue(), "")

    def test_process_sam_entry_with_failing_seq(self):
        self._prepare_for_processing_test()
        self.entry['sequence'] = "AAAAAAAATT"
        self.sam_nucl_perc_filter._process_sam_entry(self.entry)
        self.assertEqual(self.sam_nucl_perc_filter.counter_sufficient, 0)
        self.assertEqual(self.sam_nucl_perc_filter.counter_insufficient, 1)
        self.assertTrue(
            self.sam_nucl_perc_filter.output_fh_insufficient.getvalue().startswith(
                "test_seq"))
        self.assertEqual(
            self.sam_nucl_perc_filter.output_fh_sufficient.getvalue(), "")

    def test_process_sam_entry_with_passing_seq_minus_strand(self):
        self._prepare_for_processing_test()
        self.entry['sequence'] = "GGGGGGCCCCCCAAAAAAAAT"
        self.sam_nucl_perc_filter._process_sam_entry(self.entry)
        self.assertEqual(self.sam_nucl_perc_filter.counter_sufficient, 1)
        self.assertEqual(self.sam_nucl_perc_filter.counter_insufficient, 0)
        self.assertTrue(
            self.sam_nucl_perc_filter.output_fh_sufficient.getvalue().startswith(
                "test_seq"))
        self.assertEqual(
            self.sam_nucl_perc_filter.output_fh_insufficient.getvalue(), "")

    def test_process_sam_entry_with_failing_seq_minus_strand(self):
        self._prepare_for_processing_test()
        self.entry['sequence'] = "TTTTTTTTAAAAA"
        self.entry['flag'] = 16
        self.sam_nucl_perc_filter._process_sam_entry(self.entry)
        self.assertEqual(self.sam_nucl_perc_filter.counter_sufficient, 0)
        self.assertEqual(self.sam_nucl_perc_filter.counter_insufficient, 1)
        self.assertTrue(
            self.sam_nucl_perc_filter.output_fh_insufficient.getvalue().startswith(
                "test_seq"))
        self.assertEqual(
            self.sam_nucl_perc_filter.output_fh_sufficient.getvalue(), "")

    def _prepare_for_processing_test(self):
        self.entry = self._create_test_sam_enty()
        self.sam_nucl_perc_filter.counter_sufficient = 0
        self.sam_nucl_perc_filter.counter_insufficient = 0
        self.sam_nucl_perc_filter.output_fh_insufficient = StringIO()
        self.sam_nucl_perc_filter.output_fh_sufficient = StringIO()        

    def _prepare_for_processing_test(self):
        self.entry = self._create_test_sam_enty()
        self.sam_nucl_perc_filter.counter_sufficient = 0
        self.sam_nucl_perc_filter.counter_insufficient = 0
        self.sam_nucl_perc_filter.output_fh_insufficient = StringIO()
        self.sam_nucl_perc_filter.output_fh_sufficient = StringIO()        


    def _create_test_sam_enty(self):
        return({'mate_start': 0, 'mismatches': 'MD:Z:7G8', 
                'reference': 'test_org', 'sequence': 'ATCGATCGATCGATCGATCGATCG', 
                'cigar': '16M', 'flag': 0, 'query': 'test_seq', 
                'template_length': 0, 'number_of_hits': 'NH:i:1', 
                'distance': 'NM:i:1', 'phred_quality': '*', 
                'start': 4126964, 'mate': '*', 'mapping_quality': 255})

if __name__ == "__main__":
    unittest.main()

