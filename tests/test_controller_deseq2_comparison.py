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
        self.segemehl_bin = "segemehl.x"
        self.threads = 1
        self.processes = 1
        self.check_for_existing_files = False
        # DESEQ2
        self.libs = "Co_culture_replicate_1,Co_culture_replicate_2,Co_culture_replicate_3,Harboring_replicate_1,Harboring_replicate_2,Harboring_replicate_3,Infected_replicate_1,Infected_replicate_2,Infected_replicate_3,Steady_state_replicate_1,Steady_state_replicate_2,Steady_state_replicate_3,Uninfected_replicate_1,Uninfected_replicate_2,Uninfected_replicate_3"
        self.conditions = "co_culture,co_culture,co_culture,harboring,harboring,harboring,infected,infected,infected,steady_state,steady_state,steady_state,uninfected,uninfected,uninfected"
        self.replicates = "1,2,3,1,2,3,1,2,3,1,2,3,1,2,3"
        self.libs_by_species = ['human=Co_culture_replicate_1,Co_culture_replicate_2,Co_culture_replicate_3,Harboring_replicate_1,Harboring_replicate_2,Harboring_replicate_3,Infected_replicate_1,Infected_replicate_2,Infected_replicate_3,Uninfected_replicate_1,Uninfected_replicate_2,Uninfected_replicate_3', 'staphylococcus=Co_culture_replicate_1,Co_culture_replicate_2,Co_culture_replicate_3,Harboring_replicate_1,Harboring_replicate_2,Harboring_replicate_3,Infected_replicate_1,Infected_replicate_2,Infected_replicate_3,Steady_state_replicate_1,Steady_state_replicate_2,Steady_state_replicate_3']
        self.size_factor = "comparison"
        self.cooks_cutoff_off = False
        self.fc_shrinkage_off = False
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
        self.dual_fresh = "tests/test_files/reademption_analysis_dual_with_deseq_input_files_comparison_fresh"
        self.dual_fresh_output = f"{self.dual_fresh}/output"

        self.dual_fresh_copy = "a_test_project"
        self.dual_fresh_copy_output = f"{self.dual_fresh_copy}/output"
        self.dual_fresh_copy_output_human = f"{self.dual_fresh_copy_output}/human_deseq"
        #self.dual_fresh_copy_R_script_human = f"{self.dual_fresh_copy_output_human}/deseq_raw/deseq.R"
        self.dual_fresh_copy_output_staph = f"{self.dual_fresh_copy_output}/staphylococcus_deseq"

        self.dual_expected = "tests/test_files/reademption_analysis_dual_with_deseq_input_files_comparison_expected"
        self.dual_expected_output = f"{self.dual_expected}/output"
        self.dual_expected_output_human = f"{self.dual_expected_output}/human_deseq"
        #self.dual_expected_R_script_human = f"{self.dual_expected_output_human}/deseq_raw/deseq.R"
        self.dual_expected_output_staph = f"{self.dual_expected_output}/staphylococcus_deseq"

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



class TestControllerCompareWithDeseq(TestController):
    def test_deseq_comparison(self):
        self._version = 2.0
        self._copy_fresh_dual_input_files_and_structure()
        self.controller.compare_with_deseq()
        # Collect created files
        created_files = []
        for created_subfolder_files in os.walk(self.dual_fresh_copy_output_human):
            for created_subfolder_file in created_subfolder_files[2]:
                created_files.append(created_subfolder_file)

        # Collect expected files
        expected_files = []
        for expected_subfolder_files in os.walk(self.dual_expected_output_human):
            for expected_subfolder_file in expected_subfolder_files[2]:
                expected_files.append(expected_subfolder_file)
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
        #print(expected_files_and_checksums)
        #print("Created:")
        #pprint.pprint(created_files_and_checksums)
        assert expected_files == created_files

        # Staphylococcus
        # Collect the checksums of the created files
        created_files_and_checksums_staph = {}
        created_files = []
        for created_subfolder_files in os.walk(self.dual_fresh_copy_output_human):
            for created_subfolder_file in created_subfolder_files[2]:
                created_files.append(created_subfolder_file)
        # Collect the checksums of the expected files
        expected_files = []
        for expected_subfolder_files in os.walk(self.dual_expected_output_human):
            for expected_subfolder_file in expected_subfolder_files[2]:
                expected_files.append(expected_subfolder_file)
        assert expected_files == created_files

if __name__ == "__main__":
    unittest.main()
