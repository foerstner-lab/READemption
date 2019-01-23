from reademptionlib.controller_projectcreator import CreateProject
from reademptionlib.controller_alignment import PerformAlignment
from reademptionlib.controller_coverage import CalculateCoverage
from reademptionlib.controller_genequanti import GeneQuantification
from reademptionlib.controller_deseq import RunDeseq


def data_controllers():
    global test_project_name
    global version
    global project_creator
    global controller_align
    global controller_coverage
    global controller_genequanti
    global controller_deseq
    global genome_fasta
    global read_fasta_1
    global read_fasta_2
    global sam_content_1
    global sam_content_2
    global gff_content_1
    # global gff_content_2
    global overlap_output_1
    global overlap_output_2
    arg_mock_align = ArgMockAlignment()
    arg_mock_cov = ArgMockCoverage()
    arg_mock_quanti = ArgMockQuanti()
    arg_mock_deseq = ArgMockDESeq()
    version = "0.4.4.dev"
    test_project_name = "a_test_project"
    project_creator = CreateProject(arg_mock_align)
    controller_align = PerformAlignment(arg_mock_align)
    controller_coverage = CalculateCoverage(arg_mock_cov)
    controller_genequanti = GeneQuantification(arg_mock_quanti)
    controller_deseq = RunDeseq(arg_mock_deseq)
   
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
AACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGAACAGTGCGG
CGGTTGAAAATGGTTGTCGAACAAGAATTCGCTCAGATCAAACATGTTCTGCATGGTATC
ATGTCGATCGCGATTATGGCGGGACTTCTGGAGGCGCGTGGGCATCGCGTCACGGTGATC
AGGCAAGGGCAGGTAGCGACCGTACTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
TTGTCGAACAAGAATTCGCTCAGATCAAAAAAAAAAAAGGGGGTGTAAAAAAAGTGTAAA
GTGGGGTGGGTAGAGAGAGAGATTTTTTTGAGAGAGAGAAGGGTTTTTAGAGTAGAGAGG
CGCCAGCCAGATCCCGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
GGCCATTACCTTGAATCTACCGTCGATATCGCGGAATCGACTCGCCGTATCGAAAAAAAA
AAAGGGACTTCTGGAGGCGCGTGGGCATCGCGTCACGGTGAAAAAAAAAAAAAAAAAAAA
TCTGGAGGCGCGTGGGCATCGCGTCACGGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
GAATCGACTCGCCGTATCGCCGCCAGCCAGATCCCGGCCGATCAGATGATCCTGATGGCG
ATGGCGGGACTTCTGGAGGCGCGTGGGCATCGCGTCACGGTGATCAAAAAAAAAAAAAAA
GGTCAGTGCCCGGATAGCATCAACGCCGCGCTGATTTGCAAAAAAAAAAAAAAAAAAAAA
AAGTTTTTTTGTGAGAGAGAAGTTTTGAGAGAGAGTTAGAGGAAAAAAAAAAAAAAAAAA
CGCCAGCAGCACATGAACAAGTTTCGGAATGTGATCAATTTAAAAATTTATTGACTTAGG
CGCCAGCAGCACATGAACAAGTTTCGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
ATGAACAAGTTTCGGAATGTGATCAATTTAAAAATTTATTGACTTAGGAAAAAAAAAAAA
TGTGATCAATTTAAAAATTTATTGACTTAGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
GGCCATGACCTTGAATCTACCGTCGATATCGCGGAATCGACTCGCCGTATCGAAAAAAAA
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
read_01	0	SL1344	50	255	10M	*	0	0	ACAACATCCA	*	NM:i:0	MD:Z:10	NH:i:1	XA:Z:Q
"""

    sam_content_2 = """@HD	VN:1.0
@SQ	SN:SL1344	LN:960
@PG	ID:segemehl	VN:0.9.4-$Rev: 316 $ ($Date: 2011-08-18 16:37:19 +0200 (Thu, 18 Aug 2011) $)
read_01	0	SL1344	100	255	10M	*	0	0	ACAACATCCA	*	NM:i:0	MD:Z:10	NH:i:1	XA:Z:Q
read_01	0	SL1344	500	255	10M	*	0	0	ACAACATCCA	*	NM:i:0	MD:Z:10	NH:i:1	XA:Z:Q
"""

    gff_content_1 = """##gff-version 3
#!gff-spec-version 1.14
#!source-version NCBI C++ formatter 0.2
##Type DNA SL1344
SL1344	EMBL	gene	99	115	.	+	.	ID=SL1344:foo;locus_tag=SL1344_0001
SL1344	EMBL	gene	99	115	.	-	.	ID=SL1344:bar;locus_tag=SL1344_0002
SL1344	EMBL	gene	110	130	.	+	.	ID=SL1344:samba;locus_tag=SL1344_0003
SL1344	EMBL	gene	109	140	.	+	.	ID=SL1344:limbo;locus_tag=SL1344_0004
SL1344	EMBL	gene	505	550	.	-	.	ID=SL1344:rumba;locus_tag=SL1344_0005
"""

    # Currently not used
    gff_content_2 = """##gff-version 3
#!gff-spec-version 1.14
#!source-version NCBI C++ formatter 0.2
##Type DNA SL1344
SL1344	EMBL	source	1	4878012	.	+	.	organism=Salmonella enterica subsp. enterica serovar Typhimurium str. SL1344;
SL1344	EMBL	gene	169	255	.	+	.	ID=SL1344:thrL;locus_tag=SL1344_0001
SL1344	EMBL	CDS	169	252	.	+	0	ID=SL1344:thrL:unknown_transcript_1;Parent=SL1344:thrL;locus_tag=SL1344_0001;
SL1344	EMBL	start_codon	169	171	.	+	0	ID=SL1344:thrL:unknown_transcript_1;Parent=SL1344:thrL;locus_tag=SL1344_0001;
SL1344	EMBL	stop_codon	253	255	.	+	0	ID=SL1344:thrL:unknown_transcript_1;Parent=SL1344:thrL;locus_tag=SL1344_0001;
SL1344	EMBL	gene	337	2799	.	+	.	ID=SL1344:thrA;locus_tag=SL1344_0002
SL1344	EMBL	CDS	337	2796	.	+	0	ID=SL1344:thrA:unknown_transcript_1;Parent=SL1344:thrA;locus_tag=SL1344_0002;
SL1344	EMBL	start_codon	337	339	.	+	0	ID=SL1344:thrA:unknown_transcript_1;Parent=SL1344:thrA;locus_tag=SL1344_0002;
SL1344	EMBL	stop_codon	2797	2799	.	+	0	ID=SL1344:thrA:unknown_transcript_1;Parent=SL1344:thrA;locus_tag=SL1344_0002;
SL1344	EMBL	misc_feature	337	2796	.	+	.	ID=SL1344:thrA:unknown_transcript_2;Parent=SL1344:thrA;locus_tag=SL1344_0002;
SL1344	EMBL	misc_feature	337	1224	.	+	.	ID=SL1344:thrA:unknown_transcript_3;Parent=SL1344:thrA;locus_tag=SL1344_0002;
SL1344	EMBL	misc_feature	349	351	.	+	.	ID=SL1344:thrA:unknown_transcript_4;Parent=SL1344:thrA;locus_tag=SL1344_0002;
"""

    overlap_output_1 = """read_01	SL1344	100	109	+	1	SL1344	EMBL	gene	99	115	.	+	.	ID=SL1344:foo;locus_tag=SL1344_0001
read_01	SL1344	100	109	+	1	SL1344	EMBL	gene	99	115	.	-	.	ID=SL1344:bar;locus_tag=SL1344_0002
read_01	SL1344	100	109	+	1	SL1344	EMBL	gene	109	140	.	+	.	ID=SL1344:limbo;locus_tag=SL1344_0004
read_01	SL1344	500	509	+	1	SL1344	EMBL	gene	505	550	.	-	.	ID=SL1344:rumba;locus_tag=SL1344_0005
"""

    overlap_output_2 = """read_01	SL1344	1	10	+	1	no_overlap
read_01	SL1344	50	59	+	1	no_overlap
"""


class ArgMockAlignment(object):
    project_path = "a_test_project"
    min_read_length = 20
    STAR_bin = "STAR"
    indexN = 10
    threads = 1
    paired_end = False
    cutadapt = False
    segemehl = False
    processes = 1
    check_for_existing_files = False
    poly_a_clipping = True
    split = False
    realign = False
    crossalign_cleaning_str = None
    min_phred_score = None
    adapter = None
    reverse_complement = False


class ArgMockCoverage(object):
    project_path = "a_test_project"
    processes = 1
    normalize_by_uniquely = False
    non_strand_specific = False
    skip_read_count_splitting = False
    unique_only = False
    coverage_style = "global"
    clip_length = 11
    check_for_existing_files = False


class ArgMockQuanti(object):
    project_path = "a_test_project"
    min_overlap = 1
    read_region = "global"
    clip_length = 1
    paired_end = False
    no_count_split_by_alignment_no = False
    no_count_splitting_by_gene_no = False
    skip_antisense = False
    non_strand_specific = False
    processes = 1
    features = None
    allowed_features = None
    unique_only = False
    pseudocounts = False
    check_for_existing_files = False


class ArgMockDESeq(object):
    project_path = "a_test_project"
    libs = "libbar.fa,libfoo.fa"
    conditions = "condition_1,condition_2"
    cooks_cutoff_off = False
    padj_cutoff = 0.05
    shape = "circle"
    alpha = 0.5
    color_sig = "red"
    color_non_sig = "black"
    glyph_size = 8
    deseq_raw_folder = "{}/output/deseq/deseq_raw".format(project_path)
    deseq_extended_folder = "{}/output/deseq/deseq_with_annotations".format(
        project_path)
    deseq_script_path = deseq_raw_folder
    deseq_pca_heatmap_path = deseq_raw_folder
    gene_wise_quanti_combined_path = deseq_raw_folder
    deseq_tmp_session_info_script = deseq_raw_folder
    deseq_session_info = deseq_raw_folder


