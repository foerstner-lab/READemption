import hashlib
import unittest

from segemehl import Segemehl

class TestSegemehl(unittest.TestCase):
    
    fasta_file_path = "/tmp/test.fa"
    index_file_path = "/tmp/test.idx"

    def setUp(self):
        self.segemehl = Segemehl(segemehl_bin="segemehl")
        self.example_data = ExampleData()

    def test_build_index(self):
        self._create_tmp_fasta_file(self.fasta_file_path)
        self.segemehl.build_index(self.fasta_file_path, self.index_file_path)
        assert(self._sha224_of_file())

    def test_map_reads(self):
        pass

    def _create_tmp_fasta_file(self, fasta_file_path):
        fasta_fh = open(fasta_file_path, "wb")
        fasta_fh.write(self.example_data.genome_fasta_1)
        fasta_fh.close()

    def _sha224_of_file(self):
        return(hashlib.sha224(open(self.fasta_file_path).read()).hexdigest())

class ExampleData(object):

    read_fasta_1 = """>read_01
AAAAATTTTTTTAAAA
>read_02
"""

    genome_fasta_1 = """>

""" 

if __name__ == "__main__":
    unittest.main()
