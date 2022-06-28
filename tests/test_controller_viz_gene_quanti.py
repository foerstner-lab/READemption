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
        self.triple_fresh = "tests/test_files/reademption_analysis_triple_after_gene_quanti_remove_cross_aligned_reads"
        self.triple_fresh_copy = "a_test_project"
        self.triple_fresh_copy_output = f"{self.triple_fresh_copy}/output"
        self.triple_fresh_copy_output_viz_gene_quanti_human = f"{self.triple_fresh_copy_output}/human_viz_gene_quanti"
        self.triple_fresh_copy_output_viz_gene_quanti_staph = f"{self.triple_fresh_copy_output}/staphylococcus_viz_gene_quanti"
        self.triple_fresh_copy_output_viz_gene_quanti_influenza = f"{self.triple_fresh_copy_output}/influenza_viz_gene_quanti"


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
    def test_viz_gene_quanti_comparison(self):
        self._version = 2.0
        self._copy_fresh_triple_input_files_and_structure()
        self.controller.viz_gene_quanti()
        human_viz_gene_quanti_files_expected = ['rna_class_sizes.pdf', 'expression_scatter_plots.pdf']
        human_viz_gene_quanti_files_calculated = os.listdir(self.triple_fresh_copy_output_viz_gene_quanti_human)
        assert sorted(human_viz_gene_quanti_files_expected) == sorted(human_viz_gene_quanti_files_calculated)

        staphyloccocus_viz_gene_quanti_files_expected = ['rna_class_sizes.pdf', 'expression_scatter_plots.pdf']
        staphyloccocus_viz_gene_quanti_files_calculated = os.listdir(self.triple_fresh_copy_output_viz_gene_quanti_staph)
        assert sorted(staphyloccocus_viz_gene_quanti_files_expected) == sorted(staphyloccocus_viz_gene_quanti_files_calculated)

        influenza_viz_gene_quanti_files_expected = ['rna_class_sizes.pdf', 'expression_scatter_plots.pdf']
        influenza_viz_gene_quanti_files_calculated = os.listdir(self.triple_fresh_copy_output_viz_gene_quanti_influenza)
        assert sorted(influenza_viz_gene_quanti_files_expected) == sorted(influenza_viz_gene_quanti_files_calculated)

if __name__ == "__main__":
    unittest.main()
