import unittest
import sys
from io import StringIO
sys.path.append(".")
from libs.seqsizefilter import SeqSizeFilter

class TestSeqSizeFilter(unittest.TestCase):

    def setUp(self):
        self.seqsizefilter = SeqSizeFilter()
        self.example_data = ExampleData()

    def test_filter(self):
        """If the number of input file is not equal the number of
        output file an exception must be raised."""
        self.assertRaises(
            Exception, self.seqsizefilter.filter, ["a"], ["a", "b"], ["a"])

    def test_filter_entries_in_file(self):
        input_fh = StringIO(self.example_data.fasta_input)
        output_long_enough_fh = StringIO()
        output_too_short_fh = StringIO()
        self.seqsizefilter._filter_entries_in_file(
            input_fh, output_long_enough_fh, output_too_short_fh, 12)
        self.assertEqual(output_long_enough_fh.getvalue(), 
                         self.example_data.fasta_long_enough_output)

class ExampleData(object):

    fasta_input = """>test_1 12 nt
ACTGACTGACTG
>test_2 10 nt
ACTGACTGAC
>test_3  6 nt
ACTGAC
>test_4 32 nt
ACTGACTGACTGACTGACTGACTGACTGACTG
>test_5 64 nt plus line break
ACTGACTGACTGACTGACTGACTGACTGACTG
ACTGACTGACTGACTGACTGACTGACTGACTG
"""
    fasta_long_enough_output = """>test_1 12 nt
ACTGACTGACTG
>test_4 32 nt
ACTGACTGACTGACTGACTGACTGACTGACTG
>test_5 64 nt plus line break
ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG
"""

    fasta_to_short = """>test_2 10 nt
ACTGACTGAC
>test_3  6 nt
ACTGAC
"""
        
if __name__ == "__main__":
    unittest.main()
