from reademptionlib.segemehl import Segemehl


def data_segemehl():
    global fasta_file_path
    global index_file_path
    global read_fasta_file_path
    global aligning_result_path
    global unmapped_reads_path
    global segemehl
    global maxDiff
    global genome_fasta_lower
    global genome_fasta_upper
    global sam_result_aligned_1
    global sam_result_aligned_2
    global sam_result_aligned_3
    global sam_result_no_aligned
    fasta_file_path = "/tmp/test.fa"
    index_file_path = "/tmp/test.idx"
    read_fasta_file_path = "/tmp/test_reads.fa"
    aligning_result_path = "/tmp/test_aligning_results.sam"
    unmapped_reads_path = "/tmp/test_unmapped_reads.fa"
    segemehl = Segemehl(segemehl_bin="segemehl.x")
    maxDiff = None
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
    sam_result_aligned_1 = """
read_01	0	SL1344	181	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:1	XI:i:0	XA:Z:Q
"""
    sam_result_aligned_2 = """
read_03	0	SL1344	301	255	20M	*	0	0	GCTTTTTTTTCGACCAGACA	*	NM:i:1	MD:Z:18G1	NH:i:1	XI:i:0	XA:Z:Q
"""

    sam_result_aligned_3 = """
read_05	0	SL1344	301	255	20M	*	0	0	GCTTTTTTTTCGACCAGTCA	*	NM:i:2	MD:Z:17A0G1	NH:i:1	XI:i:0	XA:Z:Q
"""
    sam_result_no_aligned = """
"""
