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
        self.species = ["human=Homo sapiens", "staphylococcus=Staphylococcus aureus", "influenza=Influenza A"]
        self.threads = 1
        self.processes = 1
        self.paired_end = False
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
        self.triple_fresh = "tests/test_files/reademption_analysis_triple_after_alignment"
        self.triple_fresh_copy = "a_test_project"
        self.triple_fresh_copy_output = f"{self.triple_fresh_copy}/output"
        self.triple_fresh_copy_output_viz_align_human = f"{self.triple_fresh_copy_output}/human_viz_align"
        self.triple_fresh_copy_output_viz_align_staph = f"{self.triple_fresh_copy_output}/staphylococcus_viz_align"
        self.triple_fresh_copy_output_viz_align_influenza = f"{self.triple_fresh_copy_output}/influenza_viz_align"


    def tearDown(self):
        self._remove_triple_input_files_and_structure_copy()

    def _copy_fresh_triple_input_files_and_structure(self):
        shutil.copytree(
            self.triple_fresh,
            self.triple_fresh_copy,
        )


    def _remove_triple_input_files_and_structure_copy(self):
        if os.path.exists(self.triple_fresh_copy):
            shutil.rmtree(self.triple_fresh_copy)



class TestControllerCompareWithDeseq(TestController):
    def test_viz_align_comparison(self):
        self._version = 2.0
        self._copy_fresh_triple_input_files_and_structure()
        self.controller.viz_align()
        human_viz_align_files_expected = ['human_aligned_reads.csv', 'human_aligned_reads.pdf']
        human_viz_align_files_calculated = os.listdir(self.triple_fresh_copy_output_viz_align_human)
        assert sorted(human_viz_align_files_expected) == sorted(human_viz_align_files_calculated)

        staphylococcus_viz_align_files_expected = ['staphylococcus_aligned_reads.csv', 'staphylococcus_aligned_reads.pdf']
        staphylococcus_viz_align_files_calculated = os.listdir(self.triple_fresh_copy_output_viz_align_staph)
        assert sorted(staphylococcus_viz_align_files_expected) == sorted(staphylococcus_viz_align_files_calculated)

        influenza_viz_align_files_expected = ['influenza_aligned_reads.csv', 'influenza_aligned_reads.pdf']
        influenza_viz_align_files_calculated = os.listdir(self.triple_fresh_copy_output_viz_align_influenza)
        assert sorted(influenza_viz_align_files_expected) == sorted(influenza_viz_align_files_calculated)

if __name__ == "__main__":
    unittest.main()
