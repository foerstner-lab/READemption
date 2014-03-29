import hashlib
import os
import sys
import unittest
sys.path.append(".")
from reademptionlib.segemehl import Segemehl

class TestSegemehl(unittest.TestCase):
    """Provide general functionalities for tha actuall testing classes."""
    
    fasta_file_path = "/tmp/test.fa"
    index_file_path = "/tmp/test.idx"

    def setUp(self):
        self.segemehl = Segemehl(segemehl_bin="segemehl.x")
        self.example_data = ExampleData()
        self.maxDiff = None

    def _create_tmp_fasta_file(self, fasta_file_path, content):
        fasta_fh = open(fasta_file_path, "w")
        fasta_fh.write(content)
        fasta_fh.close()

    def _sha1_of_file(self, file_path):
        fh = open(file_path, "rb")
        content = fh.read()
        fh.close()
        return hashlib.sha1(content).hexdigest()

    def _remove_files(self, *args):
        for file_path in args:
            if os.path.exists(file_path):
                os.remove(file_path)

class TestSegemehlIndexBuilding(TestSegemehl):

    def test_build_index_lower_letters(self):
        self._create_tmp_fasta_file(
            self.fasta_file_path, self.example_data.genome_fasta_lower)
        self.segemehl.build_index([self.fasta_file_path], self.index_file_path)
        self.assertEqual(self._sha1_of_file(self.index_file_path), 
                         "78668505720e53735f807bb5485b0b38cc3cbc22")
        self._remove_files(self.fasta_file_path, self.index_file_path)

    def test_build_index_lower_letters(self):
        self._create_tmp_fasta_file(
            self.fasta_file_path, self.example_data.genome_fasta_upper)
        self.segemehl.build_index([self.fasta_file_path], self.index_file_path)
        self.assertEqual(self._sha1_of_file(self.index_file_path), 
                         "78668505720e53735f807bb5485b0b38cc3cbc22")
        self._remove_files(self.fasta_file_path, self.index_file_path)

class TestSegemehlAligning(TestSegemehl):

    read_fasta_file_path = "/tmp/test_reads.fa"
    aligning_result_path = "/tmp/test_aligning_results.sam"

    def setUp(self):
        super().setUp()
        self.large_output = LargeOutput()
        # Create an index file
        self._create_tmp_fasta_file(
            self.fasta_file_path, self.example_data.genome_fasta_upper)        
        self.segemehl.build_index([self.fasta_file_path], self.index_file_path)

    def tearDown(self):
        self._remove_files(self.fasta_file_path, self.index_file_path)

    def test_align_reads_single_read_perfect_match(self):
        """
        ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT
        ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT
        """
        read_file_content = (
            ">read_01\n" +
            "ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT" 
            + "\n")
        self.assertEqual(self._align_reads_and_give_result(read_file_content),
                        self.large_output.sam_result_aligned_1)

    def test_map_reads_single_read_not_matching(self):
        """
        ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT
        |       | |||     ||    |            |   ||  |            ||
        ATGTACCACATGAGAGAGATAGAGAGAGATTGACAACCACACACGAGAGAGAGAGAGAGT
        """
        read_file_content = (
            ">read_02\n" +
            "ATGTACCACATGAGAGAGATAGAGAGAGATTGACAACCACACACGAGAGAGAGAGAGAGT" 
            + "\n")
        self.assertEqual(self._align_reads_and_give_result(read_file_content),
                         self.large_output.sam_result_no_aligned_1)

    def test_map_reads_single_read_one_mismatch(self):
        """A 20 nt long read with 1 mismatch at 95% accu should be
        mapped.

        GCTTTTTTTTCGACCAGAGA
        |||||||||||||||||| |
        GCTTTTTTTTCGACCAGACA
        """
        read_file_content = (
            ">read_03\nGCTTTTTTTTCGACCAGACA\n")
        self.assertEqual(self._align_reads_and_give_result(read_file_content),
                        self.large_output.sam_result_aligned_2)

    def test_map_reads_single_read_two_mismatches_95(self):
        """A 20 nt long read with 2 mismatches at 95% accuracy should
        not be mapped.

        GCTTTTTTTTCGACCAGAGA
        |||||||||||||||||  |
        GCTTTTTTTTCGACCAGTCA
        """
        read_file_content = (
            ">read_04\nGCTTTTTTTTCGACCAGTCA\n")
        self.assertEqual(self._align_reads_and_give_result(read_file_content),
                         self.large_output.sam_result_no_aligned_1)

    def test_map_reads_single_read_two_mismatches_90(self):
        """A 20 nt long read with 2 mismatches at 90% accuracy should
        be mapped.

        GCTTTTTTTTCGACCAGAGA
        |||||||||||||||||  |
        GCTTTTTTTTCGACCAGTCA
        """
        read_file_content = (
            ">read_05\nGCTTTTTTTTCGACCAGTCA\n")
        self.assertEqual(
            self._align_reads_and_give_result(read_file_content, accuracy=90),
            self.large_output.sam_result_aligned_3)

    def test_map_reads_single_read_three_mismatches(self):
        """A 20 nt long read with 3 mismatches at 90% accuracy should
        not be mapped.
        GCTTTTTTTTCGACCAGAGA
        ||||| |||||||||||  |
        GCTTTATTTTCGACCAGTCA
        """

        read_file_content = (
            ">read_06\nGCTTTTTTTTCGACCAGTCA\n")
        self.assertEqual(self._align_reads_and_give_result(read_file_content),
                        self.large_output.sam_result_no_aligned_1)

    def test_map_reads_single_too_short_read(self):
        """Reads that are too short should be mapped
        """
        read_file_content = (
            ">read_07\nGCTTTTTTT\n")
        self.assertEqual(self._align_reads_and_give_result(read_file_content),
                        self.large_output.sam_result_no_aligned_1)

    def _align_reads_and_give_result(self, read_file_content, **kwargs):
        """

        - read_file_content: the content of a read file (in fasta format)
        - **kwargs: are directly given to map_reads()
        """

        self._create_tmp_fasta_file(
            self.read_fasta_file_path, read_file_content)
        self.segemehl.align_reads(self.read_fasta_file_path, self.index_file_path, 
                                [self.fasta_file_path], self.aligning_result_path,
                                **kwargs)
        result_fh = open(self.aligning_result_path)
        result = result_fh.read()
        result_fh.close()
        #self._remove_files(self.read_fasta_file_path, self.aligning_result_path)
        return result

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

class LargeOutput(object):

    sam_result_aligned_1 = """@HD	VN:1.0
@SQ	SN:SL1344	LN:960
@PG	ID:segemehl	VN:0.1.7-$Rev: 403 $ ($Date: 2013-09-12 11:46:53 +0200 (Thu, 12 Sep 2013) $)	CL:segemehl.x --query /tmp/test_reads.fa --index /tmp/test.idx --database /tmp/test.fa --outfile /tmp/test_aligning_results.sam --hitstrategy 1 --accuracy 95 --evalue 5.0 --threads 1 --silent
read_01	0	SL1344	181	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:1	XI:i:0	XA:Z:Q
"""

    sam_result_aligned_2 = """@HD	VN:1.0
@SQ	SN:SL1344	LN:960
@PG	ID:segemehl	VN:0.1.7-$Rev: 403 $ ($Date: 2013-09-12 11:46:53 +0200 (Thu, 12 Sep 2013) $)	CL:segemehl.x --query /tmp/test_reads.fa --index /tmp/test.idx --database /tmp/test.fa --outfile /tmp/test_aligning_results.sam --hitstrategy 1 --accuracy 95 --evalue 5.0 --threads 1 --silent
read_03	0	SL1344	301	255	20M	*	0	0	GCTTTTTTTTCGACCAGACA	*	NM:i:1	MD:Z:18G1	NH:i:1	XI:i:0	XA:Z:Q
"""

    sam_result_aligned_3 = """@HD	VN:1.0
@SQ	SN:SL1344	LN:960
@PG	ID:segemehl	VN:0.1.7-$Rev: 403 $ ($Date: 2013-09-12 11:46:53 +0200 (Thu, 12 Sep 2013) $)	CL:segemehl.x --query /tmp/test_reads.fa --index /tmp/test.idx --database /tmp/test.fa --outfile /tmp/test_aligning_results.sam --hitstrategy 1 --accuracy 90 --evalue 5.0 --threads 1 --silent
read_05	0	SL1344	301	255	20M	*	0	0	GCTTTTTTTTCGACCAGTCA	*	NM:i:2	MD:Z:17A0G1	NH:i:1	XI:i:0	XA:Z:Q
"""

    sam_result_no_aligned_1 = """@HD	VN:1.0
@SQ	SN:SL1344	LN:960
@PG	ID:segemehl	VN:0.1.7-$Rev: 403 $ ($Date: 2013-09-12 11:46:53 +0200 (Thu, 12 Sep 2013) $)	CL:segemehl.x --query /tmp/test_reads.fa --index /tmp/test.idx --database /tmp/test.fa --outfile /tmp/test_aligning_results.sam --hitstrategy 1 --accuracy 95 --evalue 5.0 --threads 1 --silent
"""

if __name__ == "__main__":
    unittest.main()

