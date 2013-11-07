import sys
from io import StringIO
import unittest
sys.path.append(".")
from trapllib.readalignerstats import ReadAlignerStats

class TestReadAlignerStats(unittest.TestCase):

    def setUp(self):
        self.read_aligner_stats = ReadAlignerStats()
        self.example_data = ExampleData()
    
    def test_count_fasta_entries(self):
        fasta_fh_1 = StringIO(self.example_data.fasta_seqs_1)
        self.assertEqual(
            3, 
            self.read_aligner_stats._count_fasta_entries(fasta_fh_1))
        fasta_fh_2 = StringIO(self.example_data.fasta_seqs_2)
        self.assertEqual(
            5, 
            self.read_aligner_stats._count_fasta_entries(fasta_fh_2))

def _mock_count_fasta_entries(file_path):
    return "the number of entries in %s" % file_path

class ExampleData(object):

    fasta_seqs_1 = """>test_1 a random sequence
TTTAGAAATTACACA
>test_2 another random sequence
ACGAGAAATTAAATTAAATT
>test_3 another random sequence
TAGAGACATTGGATTTTATT
"""

    fasta_seqs_2 = """>test_1 a random sequence
TTTAGAAATTACACA
>test_2 another random sequence
ACGAGAAATTAAATTAAATT
>test_3 another random sequence
TAGAGACATTGGATTTTATT
>test_4 another random sequence
TAGAGACATTGGATTTTATT
>test_5 another random sequence
TAGAGACATTGGATTTTATT
"""

    stat_file_content = """\txxx\tyyy\tzzz
Number of raw reads\t8\t1\t5
Reads long enough after clipping\t5\t23\t42
Reads too short after clipping\t2\t1\t0
Total number of mapped reads\t18\t11\t7
Total number of mappings\t22\t11\t13
Number of unmappped reads\t4\t2\t4
Number of mapped reads in genomefile1\t8\t1\t3
Number of mapped reads in genomefile2\t10\t10\t4
Number of mapping in genomefile1\t4\t2\t1
Number of mapping in genomefile2\t18\t9\t12
"""

if __name__ == "__main__":
    unittest.main()
