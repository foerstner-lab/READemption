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
        self.threads = 1
        self.processes = 1
        self.max_pvalue = 0.05
        self.items = [self.project_path, self.species]

    # make mock argument object iterable like argparse.Namespace()
    def __iter__(self):
        return (i for i in ["project_path", "species"])


class TestController(unittest.TestCase):
    def setUp(self):
        arg_mock = ArgMock()
        self.test_project_name = arg_mock.project_path
        self.controller = Controller(arg_mock)
        self.maxDiff = None
        self.dual_fresh = "tests/test_files/reademption_analysis_dual_with_deseq_input_files_comparison_expected"
        #self.dual_fresh_output = f"{self.dual_fresh}/output"

        self.dual_fresh_copy = "a_test_project"
        self.dual_fresh_copy_output = f"{self.dual_fresh_copy}/output"
        self.dual_fresh_copy_output_viz_deseq_human = f"{self.dual_fresh_copy_output}/human_viz_deseq"
        #self.dual_fresh_copy_R_script_human = f"{self.dual_fresh_copy_output_human}/deseq_raw/deseq.R"
        self.dual_fresh_copy_output_viz_deseq_staph = f"{self.dual_fresh_copy_output}/staphylococcus_viz_deseq"

        #self.dual_expected = "tests/test_files/reademption_analysis_dual_with_deseq_input_files_comparison_expected"
        #self.dual_expected_output = f"{self.dual_expected}/output"
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



class TestControllerCompareWithDeseq(TestController):
    def test_viz_deseq_comparison(self):
        self._version = 2.0
        self._copy_fresh_dual_input_files_and_structure()
        self.controller.viz_deseq()
        human_viz_deseq_files_expected = ['volcano_plots_log2_fold_change_vs_adjusted_p-value.pdf', 'MA_plots.pdf', 'volcano_plots_log2_fold_change_vs_p-value.pdf']
        human_viz_desq_files_calculated = os.listdir(self.dual_fresh_copy_output_viz_deseq_human)
        assert sorted(human_viz_deseq_files_expected) == sorted(human_viz_desq_files_calculated)
        staphyloccocus_viz_deseq_files_expected = ['volcano_plots_log2_fold_change_vs_adjusted_p-value.pdf', 'MA_plots.pdf', 'volcano_plots_log2_fold_change_vs_p-value.pdf']
        staphyloccocus_viz_desq_files_calculated = os.listdir(self.dual_fresh_copy_output_viz_deseq_staph)
        assert sorted(staphyloccocus_viz_deseq_files_expected) == sorted(staphyloccocus_viz_desq_files_calculated)

if __name__ == "__main__":
    unittest.main()
