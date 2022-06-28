import os
import unittest
import shutil
import sys
import tempfile


sys.path.append(".")
from reademptionlib.pathcreator import PathCreator


class TestPaths(unittest.TestCase):

    test_folder = "/tmp/test"
    test_config = "tmp/config.json"
    test_species = [
        "human=homo sapiens",
        "staphylococcus=Staphylococcus aureus",
        "influenza=Influenza A",
    ]
    test_files = ["foo.fa", "bar.fa"]
    test_lib_names = ["foo", "bar"]

    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        project_folders = [
            "input",
            "input/reads",
            "input/human_annotations",
            "input/human_annotations/human.gff",
            "input/human_reference_sequences",
            "input/human_reference_sequences/human.fa",
            "input/influenza_annotations",
            "input/influenza_annotations/influenza.gff",
            "input/influenza_reference_sequences",
            "input/influenza_reference_sequences/influenza.fa",
            "input/staphylococcus_annotations",
            "input/staphylococcus_annotations/staphylococcus.gff",
            "input/staphylococcus_reference_sequences",
            "input/staphylococcus_reference_sequences/staphylococcus.fa",
            "output",
            "output/align",
            "output/align/alignments",
            "output/align/index",
            "output/align/processed_reads",
            "output/align/reports_and_stats",
            "output/align/unaligned_reads",
        ]
        for folder in project_folders:
            os.mkdir(os.path.join(self.test_dir, folder))
        self.pathcreator = PathCreator(self.test_dir, self.test_species)
        self.folder_names = [
            self.pathcreator.input_folder,
            self.pathcreator.output_folder,
            self.pathcreator.align_report_folder,
            self.pathcreator.raw_stat_data_folder,
            self.pathcreator.read_fasta_folder,
            self.pathcreator.ref_seq_folders_by_species,
            self.pathcreator.annotation_folders_by_species,
            self.pathcreator.read_alignment_index_folder,
            self.pathcreator.read_alignments_folder,
            self.pathcreator.processed_reads_folder,
            self.pathcreator.unaligned_reads_folder,
            # Tests for coverage and gene quanti folders need to be adjusted
            # # after changing the subcommands in regard to dual rna seq
            # self.pathcreator.coverage_raw_folder,
            # self.pathcreator.coverage_tnoar_min_norm_folder,
            # self.pathcreator.coverage_tnoar_mil_norm_folder,
            # self.pathcreator.gene_quanti_base_folder,
            # self.pathcreator.gene_wise_quanti_combined_path,
        ]
        self.pathcreator.ref_seq_folders_by_species = {
            "human": f"{self.test_dir}/input/human_reference_sequences",
            "influenza": f"{self.test_dir}/input/influenza_reference_sequences",
            "staphylococcus": f"{self.test_dir}/input/staphylococcus_reference_sequences",
        }

        self.pathcreator.annotation_paths_by_species = {
            "human": f"{self.test_dir}/input/human_annotations",
            "influenza": f"{self.test_dir}/input/influenza_annotations",
            "staphylococcus": f"{self.test_dir}/input/staphylococcus_annotations",
        }


        self.static_files = [
            self.pathcreator.read_processing_stats_path,
            self.pathcreator.read_alignments_stats_path,
            self.pathcreator.read_file_stats,
            self.pathcreator.read_alignment_stats_table_path,
            self.pathcreator.ref_seq_file_stats,
            self.pathcreator.index_path,
        ]

    def tearDown(self):
        self._remove_folder_if_exists(self.test_folder)
        # Remove temporary directory
        shutil.rmtree(self.test_dir)

    def test_set_folder_names(self):
        for folder_name in self.folder_names:
            assert folder_name != ""
            self.assertEqual(self.folder_names.count(folder_name), 1)

    def test_set_files(self):
        self.pathcreator._set_static_files
        for file_name in self.static_files:
            assert file_name != ""
            self.assertEqual(self.static_files.count(file_name), 1)

    def test_get_sorted_folder_content(self):
        self._create_folder_with_files(self.test_folder, self.test_lib_names)
        self.assertEqual(
            self.pathcreator._get_sorted_folder_content(self.test_folder),
            sorted(self.test_lib_names),
        )
        self._remove_folder_if_exists(self.test_folder)

    def test_get_read_files(self):
        self.pathcreator.read_fasta_folder = self.test_folder
        self._create_folder_with_files(self.test_folder, self.test_lib_names)
        self.assertEqual(
            self.pathcreator.get_read_files(), sorted(self.test_lib_names)
        )
        self._remove_folder_if_exists(self.test_folder)

    def test_get_ref_seq_files(self):
        ref_seq_files = self.pathcreator.get_ref_seq_files()
        ref_seq_files_expected = [
            "human.fa",
            "influenza.fa",
            "staphylococcus.fa",
        ]
        self.assertEqual(ref_seq_files, ref_seq_files_expected)

    # first the annotation files need to be returned by species. Then this
    # test has to be adjusted
    def test_get_annotation_files(self):
        expected_annoation_files = ['human.gff', 'staphylococcus.gff', 'influenza.gff']
        self.assertEqual(
            self.pathcreator.get_annotation_files(), expected_annoation_files
        )
        self._remove_folder_if_exists(self.test_folder)

    def test_set_read_files_dep_file_lists(self):
        self.pathcreator.set_read_files_dep_file_lists_single_end(
            self.test_files, self.test_lib_names
        )
        for path_list in [
            self.pathcreator.read_paths,
            self.pathcreator.processed_read_paths,
            self.pathcreator.unaligned_reads_paths,
        ]:
            assert isinstance(path_list, list)

    def test_path_list_without_appendix(self):
        self.assertEqual(
            self.pathcreator._path_list(self.test_folder, self.test_lib_names),
            ["/tmp/test/foo", "/tmp/test/bar"],
        )

    def test_path_list_with_appendix(self):
        self.assertEqual(
            self.pathcreator._path_list(
                self.test_folder, self.test_lib_names, appendix=".boing"
            ),
            ["/tmp/test/foo.boing", "/tmp/test/bar.boing"],
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
