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
        
    def test_start_project(self):
        arg_mock = ArgMock()
        arg_mock.project_name = self.test_project_name
        self.controller.start_project(arg_mock)
        self.assertEqual(
            list(os.listdir(self.test_project_name)), 
            ['rapl.config', 'input', 'output'])
        self._remove_project_folder()

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

    def _generate_input_fasta_files(self):
        genome_fh = open("%s/%s" % (
                self.controller.paths.genome_folder, "a_genome.fa"), "w")
        read_fh_1 = open("%s/%s" % (
                self.controller.paths.read_fasta_folder, "read_foo.fa"), "w")
        read_fh_2 = open("%s/%s" % (
                self.controller.paths.read_fasta_folder, "read_bar.fa"), "w")
        genome_fh.write(self.example_data.genome_fasta)
        genome_fh.close()
        read_fh_1.write(self.example_data.read_fasta_1)
        read_fh_1.close()
        read_fh_2.write(self.example_data.read_fasta_2)
        read_fh_2.close()

    def _remove_project_folder(self):
        if os.path.exists(self.test_project_name):
            shutil.rmtree(self.test_project_name)
    
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

if __name__ == "__main__":
    unittest.main()

