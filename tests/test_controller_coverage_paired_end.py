import os
import sys
import unittest
import shutil
import hashlib
import pprint

sys.path.append(".")
from reademptionlib.controller import Controller


class ArgMock:
    def __init__(self):
        self.project_path = "a_test_project"
        self.species = ["human=Homo sapiens", "staphylococcus=Staphylococcus aureus"]
        self.min_read_length = 12
        self.threads = 1
        self.processes = 1
        self.check_for_existing_files = False
        # Coverage
        self.no_fragments = False
        self.no_norm_by_fragments = True
        self.paired_end = True
        self.unique_only = False
        self.normalize_by_uniquely = False
        self.count_cross_aligned_reads = False
        self.normalize_cross_aligned_reads_included = False
        self.skip_read_count_splitting = False
        self.non_strand_specific = False
        self.coverage_style = "global"
        self.clip_length = 11
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
        self.dual_fresh = "tests/test_files/reademption_analysis_dual_with_coverage_paired_end_input_files"
        self.dual_fresh_output = f"{self.dual_fresh}/output"

        self.dual_fresh_copy = "a_test_project"
        self.dual_fresh_copy_output = f"{self.dual_fresh_copy}/output"
        #self.dual_fresh_copy_output_human = f"{self.dual_fresh_copy_output}/human_deseq"


        self.dual_expected = "tests/test_files/reademption_analysis_dual_with_coverage_paired_expected"
        self.dual_expected_output = f"{self.dual_expected}/output"
        #self.dual_expected_output_human = f"{self.dual_expected_output}/human_deseq"
        #self.dual_expected_R_script_human = f"{self.dual_expected_output_human}/deseq_raw/deseq.R"
        #self.dual_expected_output_staph = f"{self.dual_expected_output}/staphylococcus_deseq"

    def tearDown(self):
        self._remove_dual_input_files_and_structure_copy()

    def _copy_fresh_dual_input_files_and_structure(self):
        shutil.copytree(
            self.dual_fresh,
            self.dual_fresh_copy,
        )


    def _remove_dual_input_files_and_structure_copy(self):
        if os.path.exists(self.dual_fresh_copy):
            shutil.rmtree(self.dual_fresh_copy)



class TestControllerCoveragePairedEnd(TestController):
    def test_coverage_paired_end(self):
        self._version = 2.0
        self._copy_fresh_dual_input_files_and_structure()
        self.controller.create_coverage_files()
        # Collect the checksums of the created files
        created_files_and_checksums_human = {}
        for created_subfolder_files in os.walk(self.dual_fresh_copy_output):
            for created_subfolder_file in created_subfolder_files[2]:
                output_file_path = created_subfolder_files[0] + "/" + created_subfolder_file
                # exclude R session info file and heatmap since metadata of these files
                # changes the hashsums when creating a new file due to testing
                if "align" in output_file_path:
                    continue
                # exclude gitkeeps
                if ".gitkeep" in output_file_path:
                    continue
                else:
                    hash = hashlib.md5(open(output_file_path,'rb').read()).hexdigest()
                    created_files_and_checksums_human[created_subfolder_file] = hash
        # Collect the checksums of the expected files
        expected_files_and_checksums_human = {}
        for expected_subfolder_files in os.walk(self.dual_expected_output):
            for expected_subfolder_file in expected_subfolder_files[2]:
                output_file_path = expected_subfolder_files[0] + "/" + expected_subfolder_file
                # exclude R session info file and heatmap since metadata of these files
                # changes the hashsums when creating a new file due to testing
                if "align" in output_file_path:
                    continue
                # exclude gitkeeps
                if ".gitkeep" in output_file_path:
                    continue
                else:
                    hash = hashlib.md5(open(output_file_path,'rb').read()).hexdigest()
                    expected_files_and_checksums_human[expected_subfolder_file] = hash
        # Compare the checksums
        # Can be used to compare the content of twho files line wise:
        #with open(self.dual_expected_R_script_human, 'r') as expected_R_script, open(self.dual_fresh_copy_R_script_human, 'r') as calculated_R_script:
        #    expected_contend = expected_R_script.readlines()
        #    calculated_contend = calculated_R_script.readlines()
        #    for l1, l2 in zip(expected_contend, calculated_contend):
        #        #print(l1)
        #        #print(l2)
        #        if l1 != l2:
        #            print(f"expected line: '{l1}' does not match calculated line '{l2}'")
        #    assert expected_contend == calculated_contend
        #print("Expected:")
        #pprint.pprint(expected_files_and_checksums_human)
        #print("Created:")
        #pprint.pprint(created_files_and_checksums_human)

        assert expected_files_and_checksums_human == created_files_and_checksums_human

if __name__ == "__main__":
    unittest.main()
