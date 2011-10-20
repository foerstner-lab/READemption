import os
import sys
import unittest
import shutil
sys.path.append(".")
from libs.controller import Controller

class ArgMock(object):
    project_name = None

class TestController(unittest.TestCase):

    def setUp(self):
        self.controller = Controller()
        self.test_project_name = "a_test_project"
        self.example_data = ExampleData()

    def tearDown(self):
        self._remove_project_folder()

    def _generate_input_fasta_files(self):
        genome_fh = open("%s/%s" % (
                self.controller.paths.genome_folder, "agenome.fa"), "w")
        read_fh_1 = open("%s/%s" % (
                self.controller.paths.read_fasta_folder, "libfoo.fa"), "w")
        read_fh_2 = open("%s/%s" % (
                self.controller.paths.read_fasta_folder, "libbar.fa"), "w")
        genome_fh.write(self.example_data.genome_fasta)
        genome_fh.close()
        read_fh_1.write(self.example_data.read_fasta_1)
        read_fh_1.close()
        read_fh_2.write(self.example_data.read_fasta_2)
        read_fh_2.close()

    def _generate_mapping_files(self):
        self.controller.paths.set_read_files_dep_file_lists(
            ["libfoo.fa", "libbar.fa"], 12)
        
        for file_path, sam_content in zip(
            self.controller.paths.read_mapping_result_paths, 
            [self.example_data.sam_content_1, 
             self.example_data.sam_content_2]):
            mapping_fh = open(file_path, "w")
            mapping_fh.write(sam_content)
            mapping_fh.close()
        
    def _remove_project_folder(self):
        if os.path.exists(self.test_project_name):
            shutil.rmtree(self.test_project_name)

@unittest.skip("tmp")        
class TestControllerStartProject(TestController):

    def test_start_project(self):
        arg_mock = ArgMock()
        arg_mock.project_name = self.test_project_name
        self.controller.start_project(arg_mock)
        self.assertEqual(
            list(os.listdir(self.test_project_name)), 
            ['rapl.config', 'input', 'output'])
        self._remove_project_folder()

@unittest.skip("tmp")        
class TestControllerReadMapping(TestController):

    def test_read_mapping(self):
        arg_mock = ArgMock()
        arg_mock.project_name = self.test_project_name
        self.controller.start_project(arg_mock)
        self.controller.paths._set_folder_names(
            base_path=self.test_project_name)
        self.controller.paths._set_static_file_names()
        self._generate_input_fasta_files()
        # If number of reads is less than the number of threads
        # segemehl stops. So set the number of threads to 1
        self.controller.parameters.segemehl_number_of_threads = 1
        self.controller.map_reads()
        self._remove_project_folder()

class TestControllerGRCreation(TestController):

    def test_create_gr_files(self):
        arg_mock = ArgMock()
        arg_mock.project_name = self.test_project_name
        self.controller.start_project(arg_mock)
        self.controller.paths._set_folder_names(
            base_path=self.test_project_name)
        self.controller.paths._set_static_file_names()
        self._generate_input_fasta_files()
        self._generate_mapping_files()
        self.controller.create_gr_files()
        self.assertEqual(
            sorted(list(os.listdir(self.controller.paths.gr_folder))), 
            ["libbar.fa_in_agenome.fa.minus_strand.gr",
             "libbar.fa_in_agenome.fa.plus_strand.gr",
             "libfoo.fa_in_agenome.fa.minus_strand.gr",
             "libfoo.fa_in_agenome.fa.plus_strand.gr"])
        self.assertEqual(
            sorted(list(os.listdir(
                        self.controller.paths.gr_folder_read_normalized))),
            ["libbar.fa_in_agenome.fa_norm_by_2.0_mult_by_2.0.minus_strand.gr",
             "libbar.fa_in_agenome.fa_norm_by_2.0_mult_by_2.0.plus_strand.gr",
             "libfoo.fa_in_agenome.fa_norm_by_2.0_mult_by_2.0.minus_strand.gr",
             "libfoo.fa_in_agenome.fa_norm_by_2.0_mult_by_2.0.plus_strand.gr"])
        self._remove_project_folder()
    
class ExampleData(object):

    genome_fasta = """>SL1344 genome sequence
AGAGATTACGTCTGGTTGCAAGAGATCATGACAGGGGGAATTGGTTGAAAATAAATATAT
CGCCAGCAGCACATGAACAAGTTTCGGAATGTGATCAATTTAAAAATTTATTGACTTAGG
CGGGCAGATACTTTAACCAATATAGGAATACAAGACAGACAAATAAAAATGACAGAGTAC
ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT
AACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGAACAGTGCGG
GCTTTTTTTTCGACCAGAGATCACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGT
ACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCC
AGGCAAGGGCAGGTAGCGACCGTACTTTCCGCCCCCGCGAAAATTACCAACCATCTGGTG
GCAATGATTGAAAAAACTATCGGCGGCCAGGATGCTTTGCCGAATATCAGCGATGCAGAA
CGTATTTTTTCTGACCTGCTCGCAGGACTTGCCAGCGCGCAGCCGGGATTCCCGCTTGCA
CGGTTGAAAATGGTTGTCGAACAAGAATTCGCTCAGATCAAACATGTTCTGCATGGTATC
AGCCTGCTGGGTCAGTGCCCGGATAGCATCAACGCCGCGCTGATTTGCCGTGGCGAAAAA
ATGTCGATCGCGATTATGGCGGGACTTCTGGAGGCGCGTGGGCATCGCGTCACGGTGATC
GATCCGGTAGAAAAATTGCTGGCGGTGGGCCATTACCTTGAATCTACCGTCGATATCGCG
GAATCGACTCGCCGTATCGCCGCCAGCCAGATCCCGGCCGATCACATGATCCTGATGGCG
GGCTTTACCGCCGGTAATGAAAAGGGTGAACTGGTGGTGCTGGGCCGTAATGGTTCCGAC
""" 

    read_fasta_1 = """>read_01
AACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGAACAGTGCGG
>read_02
CGGTTGAAAATGGTTGTCGAACAAGAATTCGCTCAGATCAAACATGTTCTGCATGGTATC
>read_03
ATGTCGATCGCGATTATGGCGGGACTTCTGGAGGCGCGTGGGCATCGCGTCACGGTGATC
>read_04
AGGCAAGGGCAGGTAGCGACCGTACTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>read_05
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>read_06
TTGTCGAACAAGAATTCGCTCAGATCAAAAAAAAAAAAGGGGGTGTAAAAAAAGTGTAAA
>read_07
GTGGGGTGGGTAGAGAGAGAGATTTTTTTGAGAGAGAGAAGGGTTTTTAGAGTAGAGAGG
>read_08
CGCCAGCCAGATCCCGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>read_09
GGCCATTACCTTGAATCTACCGTCGATATCGCGGAATCGACTCGCCGTATCGAAAAAAAA
>read_10
AAAGGGACTTCTGGAGGCGCGTGGGCATCGCGTCACGGTGAAAAAAAAAAAAAAAAAAAA
"""

    read_fasta_2 = """>read_01
TCTGGAGGCGCGTGGGCATCGCGTCACGGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>read_02
GAATCGACTCGCCGTATCGCCGCCAGCCAGATCCCGGCCGATCAGATGATCCTGATGGCG
>read_03
ATGGCGGGACTTCTGGAGGCGCGTGGGCATCGCGTCACGGTGATCAAAAAAAAAAAAAAA
>read_04
GGTCAGTGCCCGGATAGCATCAACGCCGCGCTGATTTGCAAAAAAAAAAAAAAAAAAAAA
>read_05
AAGTTTTTTTGTGAGAGAGAAGTTTTGAGAGAGAGTTAGAGGAAAAAAAAAAAAAAAAAA
>read_06
CGCCAGCAGCACATGAACAAGTTTCGGAATGTGATCAATTTAAAAATTTATTGACTTAGG
>read_07
CGCCAGCAGCACATGAACAAGTTTCGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>read_08
ATGAACAAGTTTCGGAATGTGATCAATTTAAAAATTTATTGACTTAGGAAAAAAAAAAAA
>read_09
TGTGATCAATTTAAAAATTTATTGACTTAGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>read_10
GGCCATGACCTTGAATCTACCGTCGATATCGCGGAATCGACTCGCCGTATCGAAAAAAAA
"""

    sam_content_1 = """@HD	VN:1.0
@SQ	SN:SL1344	LN:960
@PG	ID:segemehl	VN:0.9.4-$Rev: 316 $ ($Date: 2011-08-18 16:37:19 +0200 (Thu, 18 Aug 2011) $)
read_01	0	SL1344	1	255	10M	*	0	0	ACAACATCCA	*	NM:i:0	MD:Z:10	NH:i:1	XA:Z:Q
read_01	0	SL1344	1	255	10M	*	0	0	ACAACATCCA	*	NM:i:0	MD:Z:10	NH:i:1	XA:Z:Q
"""

    sam_content_2 = """@HD	VN:1.0
@SQ	SN:SL1344	LN:960
@PG	ID:segemehl	VN:0.9.4-$Rev: 316 $ ($Date: 2011-08-18 16:37:19 +0200 (Thu, 18 Aug 2011) $)
read_01	0	SL1344	20	255	10M	*	0	0	ACAACATCCA	*	NM:i:0	MD:Z:10	NH:i:1	XA:Z:Q
read_01	0	SL1344	20	255	10M	*	0	0	ACAACATCCA	*	NM:i:0	MD:Z:10	NH:i:1	XA:Z:Q
"""

if __name__ == "__main__":
    unittest.main()

