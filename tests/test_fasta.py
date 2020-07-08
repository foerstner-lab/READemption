from io import StringIO
import unittest
import sys

sys.path.append(".")
from reademptionlib.fasta import FastaParser


class TestFastaParser(unittest.TestCase):
    def setUp(self):
        self.fasta_parser = FastaParser()
        self.example_data = ExampleData()

    def test_parse_1(self):
        fasta_fh = StringIO(self.example_data.fasta_seqs_1)
        self.assertEqual(
            list(self.fasta_parser.entries(fasta_fh)),
            [
                ("test_1 a random sequence", "TTTAGAAATTACACA"),
                ("test_2 another random sequence", "ACGAGAAATTAAATTAAATT"),
                ("test_3 another random sequence", "TAGAGACATTGGATTTTATT"),
            ],
        )

    def test_parse_empty_file(self):
        fasta_fh = StringIO("")
        self.assertEqual(list(self.fasta_parser.entries(fasta_fh)), [])

    def test_single_entry_file_header(self):
        fasta_fh = StringIO(self.example_data.fasta_seqs_2)
        self.assertEqual(
            self.fasta_parser.single_entry_file_header(fasta_fh),
            "test_4 a random sequence",
        )

    def test_header_id_1(self):
        self.assertEqual(
            self.fasta_parser.header_id("seq_10101 An important protein"),
            "seq_10101",
        )

    def test_header_id_2(self):
        self.assertEqual(
            self.fasta_parser.header_id("seq_10101\tAn important protein"),
            "seq_10101",
        )


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

    fasta_seqs_2 = """>test_4 a random sequence
TTTAG
AAATT
ACACA
"""


if __name__ == "__main__":
    unittest.main()
