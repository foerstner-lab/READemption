import hashlib
import os
import sys
import unittest
import pysam

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
            self.fasta_file_path, self.example_data.genome_fasta_lower
        )
        self.segemehl.build_index([self.fasta_file_path], self.index_file_path)
        self.assertEqual(
            self._sha1_of_file(self.index_file_path),
            "78668505720e53735f807bb5485b0b38cc3cbc22",
        )
        self._remove_files(self.fasta_file_path, self.index_file_path)

    def test_build_index_lower_letters(self):
        self._create_tmp_fasta_file(
            self.fasta_file_path, self.example_data.genome_fasta_upper
        )
        self.segemehl.build_index([self.fasta_file_path], self.index_file_path)
        self.assertEqual(
            self._sha1_of_file(self.index_file_path),
            "78668505720e53735f807bb5485b0b38cc3cbc22",
        )
        self._remove_files(self.fasta_file_path, self.index_file_path)


class TestSegemehlAligning(TestSegemehl):

    read_fasta_file_path = "/tmp/test_reads.fa"
    aligning_result_path_bam = "/tmp/test_aligning_results.bam"
    aligning_result_path_sam = "/tmp/test_aligning_results.sam"

    def setUp(self):
        super().setUp()
        self.large_output = LargeOutput()
        # Create an index file
        self._create_tmp_fasta_file(
            self.fasta_file_path, self.example_data.genome_fasta_upper
        )
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
            ">read_01\n"
            + "ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT"
            + "\n"
        )
        self.assertEqual(
            self._remove_header_lines(self._align_reads_and_give_result(read_file_content)),
            self._remove_header_lines(self.large_output.sam_result_aligned_1),
        )

    def test_map_reads_single_read_not_matching(self):
        """
        ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT
        |       | |||     ||    |            |   ||  |            ||
        ATGTACCACATGAGAGAGATAGAGAGAGATTGACAACCACACACGAGAGAGAGAGAGAGT
        """
        read_file_content = (
            ">read_02\n"
            + "ATGTACCACATGAGAGAGATAGAGAGAGATTGACAACCACACACGAGAGAGAGAGAGAGT"
            + "\n"
        )
        self.assertEqual(
            self._remove_header_lines(self._align_reads_and_give_result(read_file_content)),
            self._remove_header_lines(self.large_output.sam_result_no_aligned_1),
        )

    def test_map_reads_single_read_one_mismatch(self):
        """A 20 nt long read with 1 mismatch at 95% accuracy should be
        mapped.

        GCTTTTTTTTCGACCAGAGA
        |||||||||||||||||| |
        GCTTTTTTTTCGACCAGACA
        """
        read_file_content = ">read_03\nGCTTTTTTTTCGACCAGACA\n"
        self.assertEqual(
            self._remove_header_lines(self._align_reads_and_give_result(read_file_content)),
            self._remove_header_lines(self.large_output.sam_result_aligned_2),
        )

    def test_map_reads_single_read_two_mismatches_95(self):
        """A 20 nt long read with 2 mismatches at 95% accuracy should
        not be mapped.

        GCTTTTTTTTCGACCAGAGA
        |||||||||||||||||  |
        GCTTTTTTTTCGACCAGTCA
        """
        read_file_content = ">read_04\nGCTTTTTTTTCGACCAGTCA\n"
        self.assertEqual(
            self._remove_header_lines(self._align_reads_and_give_result(read_file_content)),
            self._remove_header_lines(self.large_output.sam_result_no_aligned_1),
        )

    def test_map_reads_single_read_two_mismatches_90(self):
        """A 20 nt long read with 2 mismatches at 90% accuracy should
        be mapped.

        GCTTTTTTTTCGACCAGAGA
        |||||||||||||||||  |
        GCTTTTTTTTCGACCAGTCA
        """
        read_file_content = ">read_05\nGCTTTTTTTTCGACCAGTCA\n"
        print(self._remove_header_lines(self._align_reads_and_give_result(read_file_content, accuracy=90)))
        print(self._remove_header_lines(self._align_reads_and_give_result(read_file_content, accuracy=90)))
        self.assertEqual(
            self._remove_header_lines(self._align_reads_and_give_result(read_file_content, accuracy=90)),
            self._remove_header_lines(self._align_reads_and_give_result(read_file_content, accuracy=90)),
        )

    def test_map_reads_single_read_three_mismatches(self):
        """A 20 nt long read with 3 mismatches at 90% accuracy should
        not be mapped.
        GCTTTTTTTTCGACCAGAGA
        ||||| |||||||||||  |
        GCTTTATTTTCGACCAGTCA
        """

        read_file_content = ">read_06\nGCTTTTTTTTCGACCAGTCA\n"
        self.assertEqual(
            self._remove_header_lines(self._align_reads_and_give_result(read_file_content)),
            self._remove_header_lines(self.large_output.sam_result_no_aligned_1),
        )

    def test_map_reads_single_too_short_read(self):
        """Reads that are too short should be mapped
        """
        read_file_content = ">read_07\nGCTTTTTTT\n"
        self.assertEqual(
            self._remove_header_lines(self._align_reads_and_give_result(read_file_content)),
            self._remove_header_lines(self.large_output.sam_result_no_aligned_1),
        )

    def _align_reads_and_give_result(self, read_file_content, **kwargs):
        """

        - read_file_content: the content of a read file (in fasta format)
        - **kwargs: are directly given to map_reads()
        """

        self._create_tmp_fasta_file(
            self.read_fasta_file_path, read_file_content
        )
        self.segemehl.align_reads(
            self.read_fasta_file_path,
            self.index_file_path,
            [self.fasta_file_path],
            self.aligning_result_path_bam,
            **kwargs
        )
        self._convert_bam_to_sam()
        result_fh = open(self.aligning_result_path_sam)
        result = result_fh.read()
        result_fh.close()
        # self._remove_files(self.read_fasta_file_path, self.aligning_result_path)
        return result

    def _convert_bam_to_sam(self):
        pysam.view(
            "-h",
            "-o",
            self.aligning_result_path_sam,
            self.aligning_result_path_bam,
            catch_stdout=False,
        )

    def _remove_header_lines(self, alignment):
        alignment_without_header = ""
        for line in alignment.splitlines():
            if line.startswith("@"):
                continue
            else:
                alignment_without_header += line + "\n"
        return alignment_without_header

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

    sam_result_aligned_1 = """@HD	VN:1.0	SO:coordinate
@SQ	SN:SL1344	LN:960
@RG	ID:A1	SM:sample1	LB:library1	PU:unit1	PL:illumina
@PG	ID:segemehl	VN:0.3.4	CL:segemehl.x --query /tmp/test_reads.fa --index /tmp/test.idx --database /tmp/test.fa --outfile /tmp/test_aligning_results.bam --bamabafixoida --hitstrategy 1 --accuracy 95 --evalue 5.0 --threads 1
read_01	0	SL1344	181	0	60=	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:60	RG:Z:A1	YZ:Z:0
"""
    sam_result_aligned_2 = """@HD	VN:1.0	SO:coordinate
@SQ	SN:SL1344	LN:960
@RG	ID:A1	SM:sample1	LB:library1	PU:unit1	PL:illumina
@PG	ID:segemehl	VN:0.3.4	CL:segemehl.x --query /tmp/test_reads.fa --index /tmp/test.idx --database /tmp/test.fa --outfile /tmp/test_aligning_results.bam --bamabafixoida --hitstrategy 1 --accuracy 95 --evalue 5.0 --threads 1
read_03	0	SL1344	301	0	18=1X1=	*	0	0	GCTTTTTTTTCGACCAGACA	*	HI:i:0	NH:i:1	NM:i:1	MD:Z:18G1	RG:Z:A1	YZ:Z:0
"""

    sam_result_aligned_3 = """@HD	VN:1.0	SO:coordinate
@SQ	SN:SL1344	LN:960
@RG	ID:A1	SM:sample1	LB:library1	PU:unit1	PL:illumina
@PG	ID:segemehl	VN:0.3.4	CL:segemehl.x --query /tmp/test_reads.fa --index /tmp/test.idx --database /tmp/test.fa --outfile /tmp/test_aligning_results.bam --bamabafixoida --hitstrategy 1 --accuracy 90 --evalue 5.0 --threads 1
read_05	0	SL1344	301	1	17=2X1=	*	0	0	GCTTTTTTTTCGACCAGTCA	*	HI:i:0	NH:i:1	NM:i:2	MD:Z:17A0G1	RG:Z:A1	YZ:Z:0
"""

    sam_result_no_aligned_1 = """@HD	VN:1.0	SO:coordinate
@SQ	SN:SL1344	LN:960
@RG	ID:A1	SM:sample1	LB:library1	PU:unit1	PL:illumina
@PG	ID:segemehl	VN:0.3.4	CL:segemehl.x --query /tmp/test_reads.fa --index /tmp/test.idx --database /tmp/test.fa --outfile /tmp/test_aligning_results.bam --bamabafixoida --hitstrategy 1 --accuracy 95 --evalue 5.0 --threads 1
"""


if __name__ == "__main__":
    unittest.main()
