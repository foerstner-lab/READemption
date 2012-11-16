from io import StringIO
import unittest
import sys
sys.path.append(".")
from libs.annotationoverview import AnnotationOverview, AnnotationOverlap, AnnotationOverlapParser

class TestAnnotationOverview(unittest.TestCase):

    def setUp(self):
        self.example_data = ExampleData()
        self.annotation_overview = AnnotationOverview()

    @unittest.skip("Todo")
    def test_read_annotation_overlap_files(self):
        pass
    
    def test_read_annotation_overlap_file(self):
        pass

    def test_count_overlaps(self):
        pass

    def test_value_to_add(self):
        pass

    def test_count_overlaps_per_mapping(self):
        pass

    def test_mapping_key_string(self):
        pass

    def test_gene_key_string(self):
        pass

    def test_print_result_table(self):
        pass

    def test_print_result_table(self):
        pass

    def test_gff3_entry_gene_key_string(self):
        pass

    def test_write_header(self):
        pass

    def test_entries(self):
        pass

    def test_dict_to_entry(self):
        pass

class TestAnnotationOverlapParser(unittest.TestCase):
    
    def setUp(self):
        self.example_data = ExampleData()
        self.annotation_overlap_parser = AnnotationOverlapParser()
        
    def test_entries(self):
        overlap_fh = StringIO(self.example_data.overlap_output_1)
        overlaps = list(self.annotation_overlap_parser.entries(overlap_fh))
        self.assertListEqual(
            [type(overlap) for overlap in overlaps], [AnnotationOverlap] * 4)
        self.assertEqual(overlaps[0].sam_query_id, "read_01")
        self.assertEqual(overlaps[0].sam_reference, "SL1344")
        self.assertEqual(overlaps[0].sam_start, 100)
        self.assertEqual(overlaps[0].sam_end, 109)
        self.assertEqual(overlaps[0].sam_strand, "+")
        self.assertEqual(overlaps[0].sam_number_of_hits, 1)
        self.assertEqual(overlaps[0].gff_reference, "SL1344")
        self.assertEqual(overlaps[0].gff_source, "EMBL")
        self.assertEqual(overlaps[0].gff_feature, "gene")
        self.assertEqual(overlaps[0].gff_start, 99)
        self.assertEqual(overlaps[0].gff_end, 115)
        self.assertEqual(overlaps[0].gff_score, ".")
        self.assertEqual(overlaps[0].gff_strand, "+")
        self.assertEqual(overlaps[0].ggf_phase, ".")
        self.assertEqual(overlaps[0].gff_attibutes, 
                         "ID=SL1344:foo;locus_tag=SL1344_0001")
        self.assertEqual(overlaps[0].orientation, "s")
        self.assertEqual(overlaps[0].overlap_length, 10)

class TestAnnotationOverlap(unittest.TestCase):

    def setUp(self):
        self.annotation_overlap = AnnotationOverlap({
                "sam_query_id" : "Read1010", "sam_reference" : "a_genome", 
                "sam_start" : "15", "sam_end" : "200", "sam_strand" : "+", 
                "sam_number_of_hits" : "5", "gff_reference" : "a_genome", 
                "gff_source" : "gff_maker", "gff_feature" : "CDS", 
                "gff_start" : "100", "gff_end" : "300", "gff_score" : ".", 
                "gff_strand" : "+", "ggf_phase" : ".", 
                "gff_attibutes" : "locus_tag=yup"})

    def test_creation(self):
        self.assertEqual(self.annotation_overlap.sam_query_id, "Read1010")
        self.assertEqual(self.annotation_overlap.sam_reference, "a_genome")
        self.assertEqual(self.annotation_overlap.sam_start, 15)
        self.assertEqual(self.annotation_overlap.sam_end, 200)
        self.assertEqual(self.annotation_overlap.sam_strand, "+")
        self.assertEqual(self.annotation_overlap.sam_number_of_hits, 5)
        self.assertEqual(self.annotation_overlap.gff_reference, "a_genome")
        self.assertEqual(self.annotation_overlap.gff_source, "gff_maker")
        self.assertEqual(self.annotation_overlap.gff_feature, "CDS")
        self.assertEqual(self.annotation_overlap.gff_start, 100)
        self.assertEqual(self.annotation_overlap.gff_end, 300)
        self.assertEqual(self.annotation_overlap.gff_score, ".")
        self.assertEqual(self.annotation_overlap.gff_strand, "+")
        self.assertEqual(self.annotation_overlap.ggf_phase, ".")
        self.assertEqual(self.annotation_overlap.gff_attibutes, "locus_tag=yup")
        self.assertEqual(self.annotation_overlap.orientation, "s")
        self.assertEqual(self.annotation_overlap.overlap_length, 101)

    def test_orientation(self):
        self.annotation_overlap.sam_strand = "+"
        self.annotation_overlap.gff_strand, = "+"
        self.assertEqual(self.annotation_overlap._orientation(), "s")
        self.annotation_overlap.sam_strand = "-"
        self.annotation_overlap.gff_strand, = "-"
        self.assertEqual(self.annotation_overlap._orientation(), "s")
        self.annotation_overlap.sam_strand = "+"
        self.annotation_overlap.gff_strand, = "-"
        self.assertEqual(self.annotation_overlap._orientation(), "a")
        self.annotation_overlap.sam_strand = "-"
        self.annotation_overlap.gff_strand, = "+"
        self.assertEqual(self.annotation_overlap._orientation(), "a")

    def test_overlap_length_1(self):
        """Single nucleotide overlaps
        12345678901234567890
        ====     SAM
           |
           ===== GFF
        """
        self._set_positions(1, 5, 5, 10)
        self.assertEqual(self.annotation_overlap._overlap_length(), 1)

    def test_overlap_length_2(self):
        """Single nucleotide overlaps
        12345678901234567890
            ======     SAM
            |
        =====          GFF
        """
        self._set_positions(5, 10, 1, 5)
        self.assertEqual(self.annotation_overlap._overlap_length(), 1)

    def test_overlap_length_3(self):
        """Multiple nucleotides overlap
        12345678901234567890
        ==========     SAM
            ||||||
            ========== GFF
        """
        self._set_positions(5, 10, 5, 14)
        self.assertEqual(self.annotation_overlap._overlap_length(), 6)

    def test_overlap_length_4(self):
        """Multiple nucleotides overlap
        12345678901234567890
            =========== SAM
            ||||||
        ==========      GFF
        """
        self._set_positions(5, 15, 1, 10)
        self.assertEqual(self.annotation_overlap._overlap_length(), 6)

    def test_overlap_length_5(self):
        """Identical coordinates
        12345678901234567890
        ===== SAM
        |||||
        ===== GFF
        """
        self._set_positions(1, 5, 1, 5)
        self.assertEqual(self.annotation_overlap._overlap_length(), 5)

    def test_overlap_length_6(self):
        """SAM neste in GFF
        12345678901234567890
           ======     SAM
           ||||||
        ============  GFF
        """
        self._set_positions(4, 9, 1, 12)
        self.assertEqual(self.annotation_overlap._overlap_length(), 6)

    def test_overlap_length_7(self):
        """GFF nested in SAM
        12345678901234567890
        ============ SAM
           ||||||
           ======    GFF
        """
        self._set_positions(1, 12, 4, 9)
        self.assertEqual(self.annotation_overlap._overlap_length(), 6)


    def _set_positions(self, sam_start, sam_end, gff_start, gff_end):
        self.annotation_overlap.sam_start = sam_start
        self.annotation_overlap.sam_end = sam_end
        self.annotation_overlap.gff_start = gff_start
        self.annotation_overlap.gff_end = gff_end

class ExampleData(object):

    overlap_output_1 = """read_01	SL1344	100	109	+	1	SL1344	EMBL	gene	99	115	.	+	.	ID=SL1344:foo;locus_tag=SL1344_0001
read_01	SL1344	100	109	+	1	SL1344	EMBL	gene	99	115	.	-	.	ID=SL1344:bar;locus_tag=SL1344_0002
read_01	SL1344	100	109	+	1	SL1344	EMBL	gene	109	140	.	+	.	ID=SL1344:limbo;locus_tag=SL1344_0004
read_01	SL1344	500	509	+	1	SL1344	EMBL	gene	505	550	.	-	.	ID=SL1344:rumba;locus_tag=SL1344_0005
"""

if __name__ == "__main__":
    unittest.main()

