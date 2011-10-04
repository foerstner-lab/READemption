from io import StringIO
import unittest
import sys
sys.path.append(".")
from libs.sam import SamParser, SamEntry

class TestSamParser(unittest.TestCase):

    def setUp(self):
        self.sam_parser = SamParser()
        self.example_data = ExampleData()

    def test_entries_1(self):
        seq_fh = StringIO(self.example_data.sam_content_1)
        for entry in self.sam_parser.entries(seq_fh):
            self.assertTrue(isinstance(entry, SamEntry))

class TestSamEntry(unittest.TestCase):
    
    def test_sam_entry_creation_1(self):
        split_line = [
            "read_01", "0", "SL1344", "1", "255", "10M", "*", "0", "0", 
            "AACCTTGGAA",
            "*", "NM:i:0", "MD:Z:10",  "NH:i:2", "XA:Z:Q"]
        sam_entry = SamEntry(split_line)
        self.assertEqual(sam_entry.query_id, "read_01")
        self.assertEqual(sam_entry.flag, 0)
        self.assertEqual(sam_entry.reference, "SL1344")
        self.assertEqual(sam_entry.start, 1)
        self.assertEqual(sam_entry.mapping_quality, 255)
        self.assertEqual(sam_entry.cigar, "10M")
        self.assertEqual(sam_entry.mate_next_frag, "*")
        self.assertEqual(sam_entry.mate_start, 0)
        self.assertEqual(sam_entry.template_length, 0)
        self.assertEqual(
            sam_entry.sequence, "AACCTTGGAA")
        self.assertEqual(sam_entry.phred_quality, "*")
        self.assertEqual(sam_entry.distance, "NM:i:0")
        self.assertEqual(sam_entry.mismatches, "MD:Z:10")
        self.assertEqual(sam_entry.number_of_hits, "NH:i:2")
        self.assertEqual(sam_entry.xa, "XA:Z:Q")
        self.assertEqual(sam_entry.end, 10)
        self.assertEqual(sam_entry.strand, "+")
        self.assertEqual(sam_entry.number_of_hits_as_int, 2)

    def test_sam_entry_creation_2(self):
        split_line = [
            "read_02", "16", "foo|bar", "1000", "255", "40M", "*", "0", "0", 
            "ACAACATCCATGAACCGCATCAGCACCACCACCATTACCA",
            "*", "NM:i:0", "MD:Z:40",  "NH:i:5", "XA:Z:Q"]
        sam_entry = SamEntry(split_line)
        self.assertEqual(sam_entry.query_id, "read_02")
        self.assertEqual(sam_entry.flag, 16)
        self.assertEqual(sam_entry.reference, "foo|bar")
        self.assertEqual(sam_entry.start, 1000)
        self.assertEqual(sam_entry.mapping_quality, 255)
        self.assertEqual(sam_entry.cigar, "40M")
        self.assertEqual(sam_entry.mate_next_frag, "*")
        self.assertEqual(sam_entry.mate_start, 0)
        self.assertEqual(sam_entry.template_length, 0)
        self.assertEqual(
            sam_entry.sequence, "ACAACATCCATGAACCGCATCAGCACCACCACCATTACCA")
        self.assertEqual(sam_entry.phred_quality, "*")
        self.assertEqual(sam_entry.distance, "NM:i:0")
        self.assertEqual(sam_entry.mismatches, "MD:Z:40")
        self.assertEqual(sam_entry.number_of_hits, "NH:i:5")
        self.assertEqual(sam_entry.xa, "XA:Z:Q")
        self.assertEqual(sam_entry.end, 1039)
        self.assertEqual(sam_entry.strand, "-")
        self.assertEqual(sam_entry.number_of_hits_as_int, 5)

class ExampleData(object):

    sam_content_1 = """@HD	VN:1.0
@SQ	SN:SL1344	LN:960
@PG	ID:segemehl	VN:0.9.4-$Rev: 316 $ ($Date: 2011-08-18 16:37:19 +0200 (Thu, 18 Aug 2011) $)
read_01	0	SL1344	10	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:2	XA:Z:Q
read_01	0	SL1344	1000	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:2	XA:Z:Q
read_02	0	SL1344	1000	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:3	XA:Z:Q
read_02	0	SL1344	1	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:3	XA:Z:Q
read_02	0	SL1344	1000	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:3	XA:Z:Q
read_03	16	SL1344	1500	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:3	XA:Z:Q
"""

if __name__ == "__main__":
    unittest.main()
    
