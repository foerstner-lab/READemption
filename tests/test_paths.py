import os
import unittest
import shutil
import sys

sys.path.append(".")
from reademptionlib.paths import Paths


class TestPaths(unittest.TestCase):

    test_folder = "/tmp/test"
    test_files = ["foo.fa", "bar.fa"]
    test_lib_names = ["foo", "bar"]

    def setUp(self):
        self.paths = Paths(base_path=self.test_folder)
        self.folder_names = [
            self.paths.input_folder,
            self.paths.output_folder,
            self.paths.align_report_folder,
            self.paths.raw_stat_data_folder,
            self.paths.read_fasta_folder,
            self.paths.ref_seq_folder,
            self.paths.annotation_folder,
            self.paths.read_alignment_index_folder,
            self.paths.read_alignments_folder,
            self.paths.processed_reads_folder,
            self.paths.unaligned_reads_folder,
            self.paths.coverage_raw_folder,
            self.paths.coverage_tnoar_min_norm_folder,
            self.paths.coverage_tnoar_mil_norm_folder,
            self.paths.gene_quanti_base_folder,
            self.paths.gene_wise_quanti_combined_path,
        ]
        self.static_files = [
            self.paths.read_processing_stats_path,
            self.paths.read_alignments_stats_path,
            self.paths.read_file_stats,
            self.paths.read_alignment_stats_table_path,
            self.paths.ref_seq_file_stats,
            self.paths.index_path,
        ]

    def tearDown(self):
        self._remove_folder_if_exists(self.test_folder)

    def test_set_folder_names(self):
        self.paths._set_folder_names()
        for folder_name in self.folder_names:
            assert folder_name != ""
            self.assertEqual(self.folder_names.count(folder_name), 1)

    def test_set_folder_names_with_base_path(self):
        self.paths._set_folder_names()
        for folder_name in self.folder_names:
            assert folder_name != ""
            self.assertEqual(self.folder_names.count(folder_name), 1)

    def test_set_files(self):
        self.paths._set_static_files
        for file_name in self.static_files:
            assert file_name != ""
            self.assertEqual(self.static_files.count(file_name), 1)

    def test_get_sorted_folder_content(self):
        self._create_folder_with_files(self.test_folder, self.test_lib_names)
        self.assertEqual(
            self.paths._get_sorted_folder_content(self.test_folder),
            sorted(self.test_lib_names),
        )
        self._remove_folder_if_exists(self.test_folder)

    def test_required_folders(self):
        self.assertEqual(len(self.paths.required_folders()), 25)

    def test_get_read_files(self):
        self.paths.read_fasta_folder = self.test_folder
        self._create_folder_with_files(self.test_folder, self.test_lib_names)
        self.assertEqual(
            self.paths.get_read_files(), sorted(self.test_lib_names)
        )
        self._remove_folder_if_exists(self.test_folder)

    def test_get_ref_seq_files(self):
        self.paths.ref_seq_folder = self.test_folder
        self._create_folder_with_files(self.test_folder, self.test_lib_names)
        self.assertEqual(
            self.paths.get_ref_seq_files(), sorted(self.test_lib_names)
        )
        self._remove_folder_if_exists(self.test_folder)

    def test_get_annotation_files(self):
        self.paths.annotation_folder = self.test_folder
        self._create_folder_with_files(self.test_folder, self.test_lib_names)
        self.assertEqual(
            self.paths.get_annotation_files(), sorted(self.test_lib_names)
        )
        self._remove_folder_if_exists(self.test_folder)

    def test_set_read_files_dep_file_lists(self):
        self.paths.set_read_files_dep_file_lists_single_end(
            self.test_files, self.test_lib_names
        )
        for path_list in [
            self.paths.read_paths,
            self.paths.processed_read_paths,
            self.paths.primary_read_aligner_bam_paths,
            self.paths.unaligned_reads_paths,
        ]:
            assert isinstance(path_list, list)

    def test_path_list_without_appendix(self):
        self.assertEqual(
            self.paths._path_list(self.test_folder, self.test_lib_names),
            ["/tmp/test/foo", "/tmp/test/bar"],
        )

    def test_path_list_with_appendix(self):
        self.assertEqual(
            self.paths._path_list(
                self.test_folder, self.test_lib_names, appendix=".boing"
            ),
            ["/tmp/test/foo.boing", "/tmp/test/bar.boing"],
        )

    def test_set_ref_seq_paths(self):
        self.paths.set_ref_seq_paths(self.test_files)
        self.assertEqual(
            self.paths.ref_seq_paths,
            [
                "%s/%s" % (self.paths.ref_seq_folder, file)
                for file in self.test_files
            ],
        )

    def _create_folder_with_files(self, folder, file_names):
        os.mkdir(self.test_folder)
        for file_name in ["foo", "bar"]:
            open("%s/%s" % (self.test_folder, file_name), "w").close()

    def _remove_folder_if_exists(self, folder):
        if os.path.exists(folder):
            shutil.rmtree(folder)


if __name__ == "__main__":
    unittest.main()
