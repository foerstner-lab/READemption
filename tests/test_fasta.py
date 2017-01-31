import sys
sys.path.append("./tests")
from io import StringIO
from reademptionlib.fasta import FastaParser


def test_fasta_parser():
    # Define some dummy data & parser object
    fasta_parser = FastaParser()
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

    # test fasta entry
    fasta_fh = StringIO(fasta_seqs_1)
    assert list(fasta_parser.entries(fasta_fh)) == [
        ('test_1 a random sequence', 'TTTAGAAATTACACA'),
        ('test_2 another random sequence', 'ACGAGAAATTAAATTAAATT'),
        ('test_3 another random sequence', 'TAGAGACATTGGATTTTATT')]

    # test empty fasta file
    fasta_empty_fh = StringIO("")
    assert list(fasta_parser.entries(fasta_empty_fh)) == []

    # test single entry file header
    fasta_header_fh = StringIO(fasta_seqs_2)
    assert fasta_parser.single_entry_file_header(
        fasta_header_fh) == "test_4 a random sequence"

    # test header id 1
    assert fasta_parser.header_id(
        "seq_10101 An important protein") == "seq_10101"

    # test header id 2
    assert fasta_parser.header_id(
        "seq_10101\tAn important protein") == "seq_10101"
