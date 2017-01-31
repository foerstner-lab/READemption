import sys
sys.path.append("./tests")
from io import StringIO
from reademptionlib.gff3 import Gff3Parser, Gff3Entry


def test_gff3_parser_and_entry():
    # Define example data & Gff3 Parser
    gff3_parser = Gff3Parser()
    gff_content = """##gff-version 3
foo01	ta_prog	gene	100	200	.	+	.	ID=bar01;Name=blub
foo01	a_prog	gene	400	500	.	+	.	ID=bar02;Name=limbo
foo02	a_prog	gene	400	500	.	+	.	ID=bar03
foo02	a_prog	CDS	400	500	.	+	.	ID=bar03;Name=bla
"""

    # Test 'test_entries'
    gff_fh = StringIO(gff_content)
    entry_seq_ids = [entry.seq_id for entry
                     in gff3_parser.entries(gff_fh)]
    assert ['foo01', 'foo01', 'foo02', 'foo02'] == entry_seq_ids

    # Test entry creation Nr. 1
    gff3entry_1 = Gff3Entry({
        "seq_id" : "IDfoobar2342",
        "source" : "a_prog",
        "feature" : "gene",
        "start" : "1984",
        "end" : "2001",
        "score" : "1323",
        "strand" : "+",
        "phase" : ".",
        "attributes" : "ID=foo00001;Name=foo"})
    assert gff3entry_1.seq_id == "IDfoobar2342"
    assert gff3entry_1.source == "a_prog"
    assert gff3entry_1.feature == "gene"
    assert gff3entry_1.start == 1984
    assert gff3entry_1.end == 2001
    assert gff3entry_1.score == "1323"
    assert gff3entry_1.strand == "+"
    assert gff3entry_1.phase == "."
    assert gff3entry_1.attributes == {"ID": "foo00001", "Name": "foo"}

    # Test entry creation Nr. 2
    gff3entry_2 = Gff3Entry({
        "seq_id" : "accession_424",
        "source" : "make",
        "feature" : "CDS",
        "start" : "101",
        "end" : "199",
        "score" : "0.888",
        "strand" : "-",
        "phase" : "1",
        "attributes" : "ID=blub1;Name=wow01"})
    assert gff3entry_2.seq_id == "accession_424"
    assert gff3entry_2.source == "make"
    assert gff3entry_2.feature == "CDS"
    assert gff3entry_2.start == 101
    assert gff3entry_2.end == 199
    assert gff3entry_2.score == "0.888"
    assert gff3entry_2.strand == "-"
    assert gff3entry_2.phase == "1"
    assert gff3entry_2.attributes == {"ID": "blub1", "Name": "wow01"}

    # Test entry creation Nr.3
    gff3entry_3 = Gff3Entry({
        "seq_id" : "accession_111",
        "source" : "make",
        "feature" : "CDS",
        "start" : "1000",
        "end" : "1",
        "score" : "",
        "strand" : "-",
        "phase" : "",
        "attributes" : "locus_tag=boing"})
    assert gff3entry_3.start == 1
    assert gff3entry_3.end == 1000

    # Test string creation
    gff3entry_str = Gff3Entry({
        "seq_id" : "accession_111",
        "source" : "make",
        "feature" : "CDS",
        "start" : "1000",
        "end" : "1",
        "score" : "0.5",
        "strand" : "-",
        "phase" : ".",
        "attributes" : "locus_tag=boing;note=zoong"})
    assert str(gff3entry_str == "accession_111\tmake\tCDS\t0.5\t-\t"
               "locus_tag=boing;note=zoong")
    
