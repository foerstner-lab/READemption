from io import StringIO
import unittest
import sys

sys.path.append(".")
from reademptionlib.fastq import FastqParser


class TestFastqParser(unittest.TestCase):
    def setUp(self):
        self.fastq_parser = FastqParser()
        self.example_data = ExampleData()

    def test_parse_1(self):
        fastq_fh = StringIO(self.example_data.fastq_seqs_1)
        self.assertEqual(
            list(self.fastq_parser.entries(fastq_fh)),
            [
                ("test_1 a random sequence", "TTTAGAAATTACACA"),
                ("test_2 another random sequence", "ACGAGAAATTAAATTAAATT"),
                ("test_3 another random sequence", "TAGAGACATTGGATTTTATT"),
            ],
        )

    def test_parse_empty_file(self):
        fastq_fh = StringIO("")
        self.assertEqual(list(self.fastq_parser.entries(fastq_fh)), [])

    def test_single_entry_file_header(self):
        fastq_fh = StringIO(self.example_data.fastq_seqs_2)
        self.assertEqual(
            self.fastq_parser.single_entry_file_header(fastq_fh),
            "test_4 a random sequence",
        )

    def test_header_id_1(self):
        self.assertEqual(
            self.fastq_parser.header_id("seq_10101 An important protein"),
            "seq_10101",
        )

    def test_header_id_2(self):
        self.assertEqual(
            self.fastq_parser.header_id("seq_10101\tAn important protein"),
            "seq_10101",
        )


class ExampleData(object):

    fastq_seqs_1 = """@test_1 a random sequence
TTTAGAAATTACACA
+
!!!!!!!!!!!!!!!
@test_2 another random sequence
ACGAGAAATTAAATTAAATT
+
@@@@@@@@@@@@@@@@@@@@
@test_3 another random sequence
TAGAGACATTGGATTTTATT
+
!!!!!!!!!!!!!!!!!!!!
"""

    fastq_seqs_2 = """@test_4 a random sequence
TTTAGAAATTACACA
+
!!!!!!!!!!!!!!!
"""


if __name__ == "__main__":
    unittest.main()
