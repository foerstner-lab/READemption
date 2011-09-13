from io import StringIO
import unittest
import sys
sys.path.append(".")
from libs.fasta import FastaParser


class TestFasta(unittest.TestCase):

    def setUp(self):
        self.fasta_parser = FastaParser()
        self.example_data = ExampleData()

    def test_parse_1(self):
        fasta_fh = StringIO(self.example_data.fasta_seqs_1)
        self.assertEqual(
            list(self.fasta_parser.entries(fasta_fh)), 
            [('test_1 a random sequence', 'TTTAGAAATTACACA'), 
             ('test_2 another random sequence', 'ACGAGAAATTAAATTAAATT'), 
             ('test_3 another random sequence', 'TAGAGACATTGGATTTTATT')])

    def test_parse_empty_file(self):
        fasta_fh = StringIO("")
        self.assertEqual(
            list(self.fasta_parser.entries(fasta_fh)), [])

class ExampleData(object):

    fasta_seqs_1 = """>test_1 a random sequence
TTTAG
AAATT
ACACA
>test_2 another random sequence
ACGAG
AAATT
AAATT
AAATT
>test_3 another random sequence
TAGAG
ACATT
GGATT
TTATT
"""

if __name__ == "__main__":
    unittest.main()
    
