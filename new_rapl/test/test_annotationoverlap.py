from io import StringIO
import unittest
import sys
sys.path.append(".")
from libs.annotationoverlap import AnnotationOverlap
from libs.intervaltree import IntervalTree

class Gff3EntryMoc(object):

    def __init__(self, seq_id, start, end):
        self.seq_id = seq_id
        self.start = start
        self.end = end

class TestAnnotationOverlap(unittest.TestCase):

    def setUp(self):
        self.annotations_overlap = AnnotationOverlap()

    def test_add_anno_to_genome_interval_gene(self):
        for gff_entry in [
            Gff3EntryMoc("bar", 1, 2), Gff3EntryMoc("bar", 4, 8),
            Gff3EntryMoc("foo", 16, 32), Gff3EntryMoc("foo", 16, 32)]:
            self.annotations_overlap._add_anno_to_genome_interval_gene(
                gff_entry)
        self.assertListEqual(
            sorted(self.annotations_overlap.genome_interval_gene.keys()), 
            ["bar", "foo"])
        self.assertListEqual(
            sorted(self.annotations_overlap.genome_interval_gene["bar"].keys()), 
            ["1-2", "4-8"])
        self.assertEqual(
            type(self.annotations_overlap.genome_interval_gene["bar"]["1-2"]), 
            list)
        self.assertEqual(
            type(self.annotations_overlap.genome_interval_gene["bar"]["1-2"][0]), 
            Gff3EntryMoc)
        self.assertEqual(
            type(self.annotations_overlap.genome_interval_gene[
                    "foo"]["16-32"][0]),Gff3EntryMoc)
        self.assertEqual(
            type(self.annotations_overlap.genome_interval_gene[
                    "foo"]["16-32"][1]),Gff3EntryMoc)

    def test_add_anno_to_genome_and_interval(self):
        for gff_entry in [
            Gff3EntryMoc("bar", 1, 2), Gff3EntryMoc("bar", 4, 8),
            Gff3EntryMoc("foo", 16, 32), Gff3EntryMoc("foo", 16, 32)]:
            self.annotations_overlap._add_anno_to_genome_and_interval(gff_entry)
        self.assertDictEqual(
            self.annotations_overlap.genomes_and_intervals, 
            {"bar": [(1, 2), (4, 8)], "foo": [(16, 32), (16, 32)]})
        
    def test_build_interval_trees(self):
        self.annotations_overlap.genomes_and_intervals = {
            "bar": [(1, 2), (4, 8)], "foo": [(16, 32), (16, 32)]}
        self.annotations_overlap._build_interval_trees()
        self.assertListEqual(
            sorted(self.annotations_overlap.genome_and_interval_trees.keys()),
            ["bar", "foo"])
        self.assertEqual(
            type(self.annotations_overlap.genome_and_interval_trees["bar"]), 
            IntervalTree)
        
    def test_sorted_start_end(self):
        self.assertEqual(
            self.annotations_overlap._sorted_start_end(4, 15), [4, 15])
        self.assertEqual(
            self.annotations_overlap._sorted_start_end(15, 4), [4, 15])
        self.assertEqual(
            self.annotations_overlap._sorted_start_end(1000, 50), [50, 1000])

    def test_genes_of_interval(self):
        self.annotations_overlap.genome_interval_gene = {}
        self.annotations_overlap.genome_interval_gene["genome_id_x"] = {}
        self.annotations_overlap.genome_interval_gene[
            "genome_id_x"]["1-10"] = "fake_gene_list"
        self.assertEqual(self.annotations_overlap._genes_of_interval(
                "genome_id_x", 1, 10), "fake_gene_list")

    def test_mod_seq_id(self):
        self.assertEqual(
            self.annotations_overlap._mod_seq_id("foobar"), "foobar")

    def test_mod_seq_id(self):
        self.assertEqual(self.annotations_overlap._mod_seq_id(
                "gi|15791399|ref|NC_002163.1|"), "NC_002163.1")

    def test_positions_to_key_string(self):
        self.assertEqual(
            self.annotations_overlap._positions_to_key_string(1, 4), "1-4")
        self.assertEqual(
            self.annotations_overlap._positions_to_key_string(9, 100), "9-100")

    def test_write_output_line_with_hit(self):
        output_fh = StringIO()
        sam_entry = MockSamEntry("foo_read", "bar_genome", 1, 100, "+", 1)
        gff_entry = MockGff3Entry(
            "bar_genome", "make", "gene", 50, 150, ".", "+", ".", 
            "ID=bar;locus_tag=foo")
        self.annotations_overlap._write_output_line_with_hit(
            sam_entry, gff_entry, output_fh)
        self.assertEqual(
            output_fh.getvalue(), 
            "\t".join(["foo_read", "bar_genome", "1", "100", "+", "1", 
                       "bar_genome", "make", "gene", "50", "150", ".", 
                       "+", ".", "ID=bar;locus_tag=foo"]) + "\n")
        
    def test_write_output_line_without_hit(self):
        output_fh = StringIO()
        sam_entry = MockSamEntry("foo_read", "bar_genome", 1, 100, "+", 1)
        self.annotations_overlap._write_output_line_without_hit(
            sam_entry, output_fh)
        self.assertEqual(
            output_fh.getvalue(), 
            "\t".join(["foo_read", "bar_genome", "1", "100", "+", "1", 
                       "no_overlap"])  + "\n")

    @unittest.skip("Todo")
    def test_search_overlaps(self):
        pass

class MockSamEntry(object):
    
    def __init__(self, query_id, reference, start, end, strand, 
                 number_of_hits_as_int):
        self.query_id = query_id
        self.reference = reference
        self.start = start
        self.end = end
        self.strand = strand
        self.number_of_hits_as_int = number_of_hits_as_int

class MockGff3Entry(object):
    
    def __init__(self, seq_id, source, feature, start, end, score, strand, phase,
                 attribute_string):
        self.seq_id = seq_id 
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attribute_string = attribute_string

if __name__ == "__main__":
    unittest.main()
