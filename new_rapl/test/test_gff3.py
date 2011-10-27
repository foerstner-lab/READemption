from io import StringIO
import unittest
import sys
sys.path.append(".")
from libs.gff3 import Gff3Parser, Gff3Entry

class TestGff3Parser(unittest.TestCase):

    def setUp(self):
        self.gff3parser = Gff3Parser()
        self.example_data = ExampleData()

    def test_entries(self):
        gff_fh = StringIO(self.example_data.gff_content_1)
        entry_seqids = [entry.seqid for entry 
                        in self.gff3parser.entries(gff_fh)]
        self.assertListEqual(
            ['foo01', 'foo01', 'foo02', 'foo02'], entry_seqids)
        
class TestGff3Entry(unittest.TestCase):

    def test_entry_creation_1(self):
        gff3entry = Gff3Entry({
                "seqid" : "IDfoobar2342",
                "source" : "a_prog",
                "feature" : "gene",
                "start" : "1984",
                "end" : "2001",
                "score" : "1323",
                "strand" : "+",
                "phase" : ".",
                "attributes" : "ID=foo00001;Name=foo"})
        self.assertEqual(gff3entry.seqid, "IDfoobar2342")
        self.assertEqual(gff3entry.source, "a_prog")
        self.assertEqual(gff3entry.feature, "gene")
        self.assertEqual(gff3entry.start, 1984)
        self.assertEqual(gff3entry.end, 2001)
        self.assertEqual(gff3entry.score, "1323")
        self.assertEqual(gff3entry.strand, "+")
        self.assertEqual(gff3entry.phase, ".")
        self.assertDictEqual(
            gff3entry.attributes, {"ID": "foo00001", "Name": "foo"})

    def test_entry_creation_2(self):
        gff3entry = Gff3Entry({
                "seqid" : "accession_424",
                "source" : "make",
                "feature" : "CDS",
                "start" : "101",
                "end" : "199",
                "score" : "0.888",
                "strand" : "-",
                "phase" : "1",
                "attributes" : "ID=blub1;Name=wow01"})
        self.assertEqual(gff3entry.seqid, "accession_424")
        self.assertEqual(gff3entry.source, "make")
        self.assertEqual(gff3entry.feature, "CDS")
        self.assertEqual(gff3entry.start, 101)
        self.assertEqual(gff3entry.end, 199)
        self.assertEqual(gff3entry.score, "0.888")
        self.assertEqual(gff3entry.strand, "-")
        self.assertEqual(gff3entry.phase, "1")
        self.assertDictEqual(
            gff3entry.attributes, {"ID": "blub1", "Name": "wow01"})

class ExampleData(object):
    
    gff_content_1 ="""##gff-version 3
foo01	ta_prog	gene	100	200	.	+	.	ID=bar01;Name=blub
foo01	a_prog	gene	400	500	.	+	.	ID=bar02;Name=limbo
foo02	a_prog	gene	400	500	.	+	.	ID=bar03
foo02	a_prog	CDS	400	500	.	+	.	ID=bar03;Name=bla
"""

if __name__ == "__main__":
	unittest.main()


