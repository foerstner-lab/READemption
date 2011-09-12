import os
import unittest
import shutil
import sys
sys.path.append(".")
from libs.paths import Paths

class TestPaths(unittest.TestCase):

    test_folder = "/tmp/test"

    def setUp(self):
        self.paths = Paths()
        self.folder_names = [
            self.paths.input_folder, self.paths.output_folder, 
            self.paths.rna_seq_folder, self.paths.genome_folder, 
            self.paths.input_file_stats_folder, self.paths.annotation_folder, 
            self.paths.gr_folder, self.paths.gr_folder_read_normalized, 
            self.paths.gr_folder_nucl_normalized, 
            self.paths.read_mapping_index_folder, 
            self.paths.annotation_hit_folder, 
            self.paths.no_annotation_hit_folder, 
            self.paths.annotation_hit_overview_folder, 
            self.paths.annotation_hit_overview_read_normalized_folder, 
            self.paths.annotation_hit_overview_nucl_normalized_folder, 
            self.paths.annotation_hit_overview_rpkm_normalized_folder, 
            self.paths.read_tracing_folder, 
            self.paths.report_folder, self.paths.read_mappings_folder, 
            self.paths.clipped_reads_folder, self.paths.unmapped_reads_folder]
        self.file_names = [
            self.paths.config_file, self.paths.read_file_stats, 
            self.paths.genome_file_stats, self.paths.annotation_file_stats, 
            self.paths.tracing_summary_file, self.paths.report_tex_file, 
            self.paths.lib_genome_read_mapping_summary, 
            self.paths.mapping_length_hist_r_file, 
            self.paths.mapping_length_hist_pdf_file]

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_set_folder_names(self):
        self.paths._set_folder_names()
        for folder_name in self.folder_names:
            assert(folder_name != '')
            self.assertEqual(self.folder_names.count(folder_name), 1)
            
    def test_set_file_names(self):
        self.paths._set_folder_names()
        for file_name in self.file_names:
            assert(file_name != '')
            self.assertEqual(self.file_names.count(file_name), 1)

    def test_get_sorted_folder_content(self):
        os.mkdir(self.test_folder)
        for file_name in ["foo", "bar"]:
            open("%s/%s" % (self.test_folder, file_name), "w").close()
        self.assertEqual(
            self.paths._get_sorted_folder_content(self.test_folder),
            ['bar', 'foo'])
        shutil.rmtree(self.test_folder)
            
if __name__ == "__main__":
    unittest.main()
