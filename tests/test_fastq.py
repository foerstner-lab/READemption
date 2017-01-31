import sys
sys.path.append("./tests")
from io import StringIO
from reademptionlib.fastq import FastqParser


def test_fastq_parser():
    # Define some dummy data & parser object
    fastq_parser = FastqParser()
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

    # test fastq entry
    fastq_fh = StringIO(fastq_seqs_1)
    assert list(fastq_parser.entries(fastq_fh)) == [
        ('test_1 a random sequence', 'TTTAGAAATTACACA'),
        ('test_2 another random sequence', 'ACGAGAAATTAAATTAAATT'),
        ('test_3 another random sequence', 'TAGAGACATTGGATTTTATT')]

    # test empty fastq file
    fastq_empty_fh = StringIO("")
    assert list(fastq_parser.entries(fastq_empty_fh)) == []

    # test single entry file header
    fastq_header_fh = StringIO(fastq_seqs_2)
    assert fastq_parser.single_entry_file_header(
        fastq_header_fh) == "test_4 a random sequence"

    # test header id 1
    assert fastq_parser.header_id(
        "seq_10101 An important protein") == ("seq_10101")

    # test header id 2
    assert fastq_parser.header_id(
        "seq_10101\tAn important protein") == ("seq_10101")

