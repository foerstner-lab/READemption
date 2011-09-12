import hashlib
import os
import sys
import unittest
sys.path.append(".")
from libs.segemehl import Segemehl

class TestSegemehl(unittest.TestCase):
    
    fasta_file_path = "/tmp/test.fa"
    index_file_path = "/tmp/test.idx"

    def setUp(self):
        self.segemehl = Segemehl(segemehl_bin="segemehl")
        self.example_data = ExampleData()

    def test_build_index_lower_letters(self):
        self._create_tmp_fasta_file(
            self.fasta_file_path, self.example_data.genome_fasta_lower)
        self.segemehl.build_index([self.fasta_file_path], self.index_file_path)
        self.assertEqual(self._sha1_of_file(self.index_file_path), 
                         "78668505720e53735f807bb5485b0b38cc3cbc22")
        self._remove_files()

    def test_build_index_lower_letters(self):
        self._create_tmp_fasta_file(
            self.fasta_file_path,self.example_data.genome_fasta_upper)
        self.segemehl.build_index([self.fasta_file_path], self.index_file_path)
        self.assertEqual(self._sha1_of_file(self.index_file_path), 
                         "78668505720e53735f807bb5485b0b38cc3cbc22")

    def test_map_reads(self):
        pass

    def _create_tmp_fasta_file(self, fasta_file_path, content):
        fasta_fh = open(fasta_file_path, "w")
        fasta_fh.write(content)
        fasta_fh.close()

    def _sha1_of_file(self, file_path):
        return(hashlib.sha1(open(file_path, "rb").read()).hexdigest())
    
    def _remove_files(self, file_paths):
        for file_path in file_paths:
            os.remove(file_path)

class ExampleData(object):

    genome_fasta_lower = """>SL1344 genome sequence
agagattacgtctggttgcaagagatcatgacagggggaattggttgaaaataaatatat
cgccagcagcacatgaacaagtttcggaatgtgatcaatttaaaaatttattgacttagg
cgggcagatactttaaccaatataggaatacaagacagacaaataaaaatgacagagtac
acaacatccatgaaccgcatcagcaccaccaccattaccaccatcaccattaccacaggt
aacggtgcgggctgacgcgtacaggaaacacagaaaaaagcccgcacctgaacagtgcgg
gcttttttttcgaccagagatcacgaggtaacaaccatgcgagtgttgaagttcggcggt
acatcagtggcaaatgcagaacgttttctgcgtgttgccgatattctggaaagcaatgcc
aggcaagggcaggtagcgaccgtactttccgcccccgcgaaaattaccaaccatctggtg
gcaatgattgaaaaaactatcggcggccaggatgctttgccgaatatcagcgatgcagaa
cgtattttttctgacctgctcgcaggacttgccagcgcgcagccgggattcccgcttgca
cggttgaaaatggttgtcgaacaagaattcgctcagatcaaacatgttctgcatggtatc
agcctgctgggtcagtgcccggatagcatcaacgccgcgctgatttgccgtggcgaaaaa
atgtcgatcgcgattatggcgggacttctggaggcgcgtgggcatcgcgtcacggtgatc
gatccggtagaaaaattgctggcggtgggccattaccttgaatctaccgtcgatatcgcg
gaatcgactcgccgtatcgccgccagccagatcccggccgatcacatgatcctgatggcg
ggctttaccgccggtaatgaaaagggtgaactggtggtgctgggccgtaatggttccgac
""" 

    genome_fasta_upper = """>SL1344 genome sequence
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


if __name__ == "__main__":
    unittest.main()
