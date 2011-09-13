import os
import unittest
import shutil
import sys
sys.path.append(".")
from libs.paths import Paths

class TestPaths(unittest.TestCase):

    test_folder = "/tmp/test"
    test_files = ["foo", "bar"]

    def setUp(self):
        self.paths = Paths()
        self.folder_names = [
            self.paths.input_folder, self.paths.output_folder, 
            self.paths.read_fasta_folder, self.paths.genome_folder, 
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
        self.static_files = [
            self.paths.config_file, self.paths.read_file_stats, 
            self.paths.genome_file_stats, self.paths.annotation_file_stats, 
            self.paths.tracing_summary_file, self.paths.report_tex_file, 
            self.paths.lib_genome_read_mapping_summary, 
            self.paths.mapping_length_hist_r_file, 
            self.paths.mapping_length_hist_pdf_file]

    def tearDown(self):
        self._remove_folder_if_exists(self.test_folder)

    def test_set_folder_names(self):
        self.paths._set_folder_names()
        for folder_name in self.folder_names:
            assert(folder_name != '')
            self.assertEqual(self.folder_names.count(folder_name), 1)

    def test_set_folder_names_with_base_path(self):
        self.paths._set_folder_names(base_path="/tmp")
        for folder_name in self.folder_names:
            assert(folder_name != '')
            self.assertEqual(self.folder_names.count(folder_name), 1)

    def test_set_file_names(self):
        self.paths._set_static_file_names
        for file_name in self.static_files:
            assert(file_name != '')
            self.assertEqual(self.static_files.count(file_name), 1)

    def test_get_sorted_folder_content(self):
        self._create_folder_with_files(self.test_folder, self.test_files)
        self.assertEqual(
            self.paths._get_sorted_folder_content(self.test_folder),
            sorted(self.test_files))
        self._remove_folder_if_exists(self.test_folder)

    def test_required_folders(self):
        self.assertEqual(len(self.paths.required_folders()), 21)

    def test_get_read_file_names(self):
        self.paths.read_fasta_folder = self.test_folder
        self._create_folder_with_files(self.test_folder, self.test_files)
        self.assertEqual(self.paths._get_read_file_names(),
                         sorted(self.test_files))
        self._remove_folder_if_exists(self.test_folder)

    def test_get_genome_file_names(self):
        self.paths.genome_folder = self.test_folder
        self._create_folder_with_files(self.test_folder, self.test_files)
        self.assertEqual(self.paths._get_genome_file_names(),
                         sorted(self.test_files))
        self._remove_folder_if_exists(self.test_folder)

    def test_set_read_files_dep_file_lists(self):
        self.paths.set_read_files_dep_file_lists(self.test_files, 20)
        for path_list in  [
            self.paths.read_file_paths, self.paths.clipped_read_file_paths, 
            self.paths.clipped_read_file_long_enough_paths, 
            self.paths.clipped_read_file_too_short_paths,
            self.paths.read_mapping_result_paths,
            self.paths.unmapped_reads_path]:
            assert(isinstance(path_list, list))
            
    def test_path_list_without_appendix(self):
        self.assertEqual(
            self.paths._path_list(self.test_folder, self.test_files),
            ['/tmp/test/foo', '/tmp/test/bar'])
        
    def test_path_list_with_appendix(self):
        self.assertEqual(
            self.paths._path_list(
                self.test_folder, self.test_files, appendix=".boing"),
            ['/tmp/test/foo.boing', '/tmp/test/bar.boing'])

    def test_set_genome_paths(self):
        self.paths.set_genome_paths(self.test_files)
        self.assertEqual(
            self.paths.genome_file_path, 
            ["%s/%s" % (self.paths.genome_folder, file) for 
             file in self.test_files])

    # def test_read_file_path(self):
    #     self.paths.read_fasta_folder = self.test_folder
    #     self._test_a_path(self.paths.read_file_path)

    # def test_read_file_paths(self):
    #     self.paths.read_fasta_folder = self.test_folder
    #     self._test_a_list_of_paths(self.paths.read_file_paths)
        
    # def test_clipped_read_file_path(self):
    #     self.paths.clipped_reads_folder = self.test_folder
    #     self._test_a_path(
    #         self.paths.clipped_read_file_path, appendix=".clipped.fa")

    # def test_clipped_read_file_paths(self):
    #     self.paths.clipped_reads_folder = self.test_folder
    #     self._test_a_list_of_paths(
    #         self.paths.clipped_read_file_paths, appendix=".clipped.fa")

    # def test_genome_file_path(self):
    #     self.paths.genome_folder = self.test_folder
    #     self._test_a_path(self.paths.genome_file_path)

    # def test_genome_file_paths(self):
    #     self.paths.genome_folder = self.test_folder
    #     self._test_a_list_of_paths(self.paths.genome_file_paths)
    
    # def _test_a_list_of_paths(self, path_list_returning_function, 
    #                           appendix=''):
    #     self.assertEqual(
    #         path_list_returning_function(self.test_files),
    #         ["%s/%s%s" % (self.test_folder, file_name, appendix) for file_name
    #          in self.test_files])

    # def _test_a_path(self, path_returning_function, appendix=''):
    #     self.assertEqual(
    #         path_returning_function(self.test_files[0]),
    #         "%s/%s%s" % (self.test_folder, self.test_files[0], appendix))

    def _create_folder_with_files(self, folder, file_names):
        os.mkdir(self.test_folder)
        for file_name in ["foo", "bar"]:
            open("%s/%s" % (self.test_folder, file_name), "w").close()

    def _remove_folder_if_exists(self, folder):
        if os.path.exists(folder):
            shutil.rmtree(folder)

if __name__ == "__main__":
    unittest.main()
