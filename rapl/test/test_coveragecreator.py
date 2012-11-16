import sys
from io import StringIO
import unittest
sys.path.append(".")
from libs.coveragecreator import CoverageCreator

class SamEntryMock(object):
    
    def __init__(self, start, end, strand, number_of_hits_as_int):
        self.start = start
        self.end = end
        self.strand = strand
        self.number_of_hits_as_int = number_of_hits_as_int
        self.reference = "test_ref_2"

class SamParserMock(object):

    def __init__(self, entries):
        self.entries = entries
    
    def entries_bam(self, bam_file):
        for entry in self.entries:
            yield(entry)
            
    def ref_seq_ids_and_lengths_bam(self, bam_file):
        return({"test_ref_1" : 10,
                "test_ref_2" : 10,})

class TestCoverageCreator(unittest.TestCase):

    def setUp(self):
        self.coverage_creator = CoverageCreator()


    def test_init_coverage_lists(self):
        self.coverage_creator._sam_parser = SamParserMock(
            self._generate_sam_entries())
        self.coverage_creator.init_coverage_lists("dummy_bam")
        self.assertEqual(
            self.coverage_creator.elements_and_coverages,
             {'minus': 
              {'test_ref_1': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
               'test_ref_2': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
               },
              'plus': 
              {'test_ref_1': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
               'test_ref_2': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
               }})

    def _generate_sam_entries(self):
        entries = []
        for start, end, strand, number_of_hits_as_int in [
            (1, 4, "+", 1),
            (4, 7, "+", 2),
            (1, 3, "-", 4)]:
            entries.append(SamEntryMock(
                    start, end, strand, number_of_hits_as_int))
        return(entries)

    def test_count_coverage(self):
        self.coverage_creator.elements_and_coverages = {
            'minus': {'test_ref_1': [0.0, 0.0, 0.0, 0.0, 0.0, 
                                     0.0, 0.0, 0.0, 0.0, 0.0],
             'test_ref_2': [0.0, 0.0, 0.0, 0.0, 0.0, 
                            0.0, 0.0, 0.0, 0.0, 0.0]},
            'plus': {'test_ref_1': [0.0, 0.0, 0.0, 0.0, 0.0, 
                                    0.0, 0.0, 0.0, 0.0, 0.0],
             'test_ref_2': [0.0, 0.0, 0.0, 0.0, 0.0, 
                            0.0, 0.0, 0.0, 0.0, 0.0]}}
        self.coverage_creator._sam_parser = SamParserMock(
            self._generate_sam_entries())
        self.coverage_creator.count_coverage("dummy_bam")
        self.assertEqual(
            self.coverage_creator.elements_and_coverages,
             {'minus': 
              {'test_ref_1': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
               'test_ref_2': [-0.25, -0.25, -0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
               },
              'plus': 
              {'test_ref_1': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
               'test_ref_2': [1.0, 1.0, 1.0, 1.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0]
               }})
        

if __name__ == "__main__":
    unittest.main()

