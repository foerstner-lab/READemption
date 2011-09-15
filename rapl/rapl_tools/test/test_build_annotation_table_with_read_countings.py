import unittest
import sys
sys.path.append('.')
sys.path.append('..')
from build_annotation_table_with_read_countings import AnnotationMappingTableBuilder


class TestAnnotationMappingTableBuilder(unittest.TestCase):

    def setUp(self):
        annotation_file = "an_annotation_file"
        overlap_files = ["foo", "bar"]
        min_overlap = 1
        strand_orientation = "s"
        normalization_factors = "1000:2000"
        rpkm = True
        output_file = "an_output_file"
        self.table_builder = AnnotationMappingTableBuilder(
            annotation_file, overlap_files, min_overlap, strand_orientation, 
            normalization_factors, rpkm, output_file)

    def test__calc_rpkm_working_1(self):
        self.assertEqual(
            # input: couting, normalization factor, gene length
            10000, self.table_builder._calc_rpkm(10, 5000, 200))

    def test__calc_rpkm_working_2(self):
        self.assertEqual(
            # input: couting, normalization factor, gene length
            0.4, self.table_builder._calc_rpkm(1, 5000000, 500))

    def test__calc_rpkm_norm_factor_zero(self):
        self.assertEqual(
            # input: couting, normalization factor, gene length
            self.table_builder._calc_rpkm(10, 0, 200), 0)

    def test__calc_rpkm_gene_length_zero(self):
        self.assertEqual(
            # input: couting, normalization factor, gene length
            self.table_builder._calc_rpkm(10, 0, 200), 0)


if __name__ == "__main__": 
    unittest.main()
