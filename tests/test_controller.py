import os
import sys
import unittest
import shutil
import hashlib

sys.path.append(".")
from reademptionlib.controller import Controller


class ArgMock:
    def __init__(self):
        self.project_path = "a_test_project"
        self.species = ["human=Homo sapiens", "staphylococcus=Staphylococcus aureus", "influenza=Influenza A"]
        self.min_read_length = 12
        self.segemehl_bin = "segemehl.x"
        self.threads = 1
        self.segemehl_accuracy = 100
        self.segemehl_evalue = 5.0
        self.paired_end = False
        self.processes = 1
        self.check_for_existing_files = False
        self.poly_a_clipping = False
        self.progress = False
        self.split = True
        self.realign = False
        self.crossalign_cleaning = None
        self.fastq = True
        self.min_phred_score = None
        self.adapter = None
        self.reverse_complement = False
        self.no_count_split_by_alignment_no = False
        self.no_count_splitting_by_gene_no = False
        self.remove_cross_aligned_reads = True
        self.non_strand_specific = False
        self.min_overlap = 1
        self.read_region = "global"
        self.clip_length = 11
        # Gene quantification
        self.no_count_split_by_alignment_no = False
        self.no_count_splitting_by_gene_no = False
        self.add_antisense = False
        self.antisense_only = False
        self.non_strand_specific = False
        self.allowed_features = None
        self.unique_only = False
        self.pseudocounts = False
        self.no_norm_by_fragments = False
        # Coverage
        self.normalize_by_uniquely = False
        self.skip_read_count_splitting = False
        self.coverage_style = "global"
        self.count_cross_aligned_reads = False
        self.normalize_cross_aligned_reads_included = False
        self.no_fragments = False
        self.items = [self.project_path, self.species, self.min_read_length]

    # make mock argument object iterable like argparse.Namespace()
    def __iter__(self):
        return (i for i in ["project_path", "species", "min_read_length"])


class TestController(unittest.TestCase):
    def setUp(self):
        arg_mock = ArgMock()
        self.test_project_name = arg_mock.project_path
        self.controller = Controller(arg_mock)
        self.maxDiff = None
        self.triple_fresh = "tests/test_files/reademption_analysis_triple_with_input_files"
        self.triple_fresh_copy = "a_test_project"
        self.triple_fresh_copy_output = "a_test_project/output"

        self.triple_expected_after_alignment = "tests/test_files/reademption_analysis_triple_after_alignment"
        self.triple_expected_after_alignment_output = "tests/test_files/reademption_analysis_triple_after_alignment/output"
        self.triple_expected_after_alignment_copy = "a_test_project"
        self.triple_expected_after_alignment_copy_output = "a_test_project/output"
        self.triple_expected_after_gene_quanti_output = "tests/test_files/reademption_analysis_triple_after_gene_quanti_remove_cross_aligned_reads/output"
        self.triple_expected_after_coverage_output = "tests/test_files/reademption_analysis_triple_after_coverage/output"

    def tearDown(self):
        self._remove_triple_input_files_and_structure_copy()
        self._remove_project_folder()

    def _copy_fresh_triple_input_files_and_structure(self):
        shutil.copytree(
            self.triple_fresh,
            self.triple_fresh_copy,
        )

    def _copy_aligned_files_and_structure(self):
        shutil.copytree(
            self.triple_expected_after_alignment,
            self.triple_expected_after_alignment_copy,
        )

    def _remove_triple_input_files_and_structure_copy(self):
        if os.path.exists(self.triple_fresh_copy):
            shutil.rmtree(self.triple_fresh_copy)
        if os.path.exists(self.triple_expected_after_alignment_copy):
            shutil.rmtree(self.triple_expected_after_alignment_copy)

    def _remove_project_folder(self):
        if os.path.exists(self.test_project_name):
            shutil.rmtree(self.test_project_name)


class TestControllerCreateProject(TestController):
    def test_create_project(self):
        self._version = 2.0
        self.controller.create_project(self._version)
        self.assertEqual(
            set(list(os.listdir(self.test_project_name))),
            set(["input", "output", "config.json"]),
        )

class TestControllerReadAligning(TestController):
    def test_read_aligning(self):
        self._version = 2.0
        self._copy_fresh_triple_input_files_and_structure()
        self.controller.align_reads()
        # Collect the checksums of the created files
        #created_files_and_checksums = {}
        created_files = []
        for created_subfolder_files in os.walk(self.triple_fresh_copy_output):
            for created_subfolder_file in created_subfolder_files[2]:
                output_file_path = created_subfolder_files[0] + "/" + created_subfolder_file
                # exclude gitkeeps
                if ".gitkeep" in output_file_path:
                    continue
                # exclude gz zipped files because of different checksums due to meta information
                if not output_file_path.endswith(".gz"):
                    #hash = hashlib.md5(open(output_file_path,'rb').read()).hexdigest()
                    #created_files_and_checksums[created_subfolder_file] = hash
                    created_files.append(created_subfolder_file)
        # Collect the checksums of the expected files
        #expected_files_and_checksums = {}
        expected_files = []
        for expected_subfolder_files in os.walk(self.triple_expected_after_alignment_output):
            for expected_subfolder_file in expected_subfolder_files[2]:
                output_file_path = expected_subfolder_files[0] + "/" + expected_subfolder_file
                # exclude gitkeeps
                if ".gitkeep" in output_file_path:
                    continue
                # exclude gz zipped files because of different checksums due to meta information
                if not output_file_path.endswith(".gz"):
                    #hash = hashlib.md5(open(output_file_path,'rb').read()).hexdigest()
                    #expected_files_and_checksums[expected_subfolder_file] = hash
                    expected_files.append(expected_subfolder_file)
        # Compare the checksums
        assert expected_files == created_files


class TestControllerGeneQuantification(TestController):
    def test_gene_quantification(self):
        self._version = 2.0
        self._copy_aligned_files_and_structure()
        self.controller.quantify_gene_wise()
        # Collect the checksums of the created files
        created_files_and_checksums = {}
        for created_subfolder_files in os.walk(self.triple_expected_after_alignment_copy_output):
            for created_subfolder_file in created_subfolder_files[2]:
                output_file_path = created_subfolder_files[0] + "/" + created_subfolder_file
                # exclude gz zipped files because of different checksums due to meta information
                if not output_file_path.endswith(".gz"):
                    hash = hashlib.md5(open(output_file_path,'rb').read()).hexdigest()
                    created_files_and_checksums[created_subfolder_file] = hash
        # Collect the checksums of the expected files
        expected_files_and_checksums = {}
        for expected_subfolder_files in os.walk(self.triple_expected_after_gene_quanti_output):
            for expected_subfolder_file in expected_subfolder_files[2]:
                output_file_path = expected_subfolder_files[0] + "/" + expected_subfolder_file
                # exclude gz zipped files because of different checksums due to meta information
                if not output_file_path.endswith(".gz"):
                    hash = hashlib.md5(open(output_file_path,'rb').read()).hexdigest()
                    expected_files_and_checksums[expected_subfolder_file] = hash
        # Compare the checksums
        assert expected_files_and_checksums == created_files_and_checksums

class TestControllerCoverage(TestController):
    def test_coverage(self):
        self._version = 2.0
        self._copy_aligned_files_and_structure()
        self.controller.create_coverage_files()
        # Collect the checksums of the created files
        created_files_and_checksums = {}
        for created_subfolder_files in os.walk(self.triple_expected_after_alignment_copy_output):
            for created_subfolder_file in created_subfolder_files[2]:
                output_file_path = created_subfolder_files[0] + "/" + created_subfolder_file
                # exclude gz zipped files because of different checksums due to meta information
                if not output_file_path.endswith(".gz"):
                    hash = hashlib.md5(open(output_file_path,'rb').read()).hexdigest()
                    created_files_and_checksums[created_subfolder_file] = hash
        # Collect the checksums of the expected files
        expected_files_and_checksums = {}
        for expected_subfolder_files in os.walk(self.triple_expected_after_coverage_output):
            for expected_subfolder_file in expected_subfolder_files[2]:
                output_file_path = expected_subfolder_files[0] + "/" + expected_subfolder_file
                # exclude gz zipped files because of different checksums due to meta information
                if not output_file_path.endswith(".gz"):
                    hash = hashlib.md5(open(output_file_path,'rb').read()).hexdigest()
                    expected_files_and_checksums[expected_subfolder_file] = hash
        # Compare the checksums
        print("Expected:")
        print(expected_files_and_checksums)
        print("Created:")
        print(created_files_and_checksums)
        assert expected_files_and_checksums == created_files_and_checksums


if __name__ == "__main__":
    unittest.main()
