import os
import sys
import unittest
import shutil
sys.path.append(".")
from trapllib.controller import Controller

class ArgMock(object):
    project_path = "a_test_project"
    min_read_length = 12
    segemehl_bin = "segemehl"
    threads = 1
    segemehl_accuracy = 95
    segemehl_evalue = 5.0

class TestController(unittest.TestCase):

    def setUp(self):
        arg_mock = ArgMock()
        self.test_project_name = arg_mock.project_path
        self.controller = Controller(arg_mock)
        self.example_data = ExampleData()
        self.maxDiff = None

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
            self.controller.paths.read_mapping_result_sam_paths, 
            [self.example_data.sam_content_1, 
             self.example_data.sam_content_2]):
            mapping_fh = open(file_path, "w")
            mapping_fh.write(sam_content)
            mapping_fh.close()

    def _generate_annotation_files(self):
        annotation_fh = open(
            "%s/some_annos.gff" % 
            self.controller.paths.annotation_folder, "w")
        print(self.controller.paths.annotation_folder)
        annotation_fh.write(self.example_data.gff_content_1)
        annotation_fh.close()
        
    def _remove_project_folder(self):
        if os.path.exists(self.test_project_name):
            shutil.rmtree(self.test_project_name)

class TestControllerStartProject(TestController):

    def test_start_project(self):
        self.controller.start_project()
        self.assertEqual(
            list(os.listdir(self.test_project_name)), 
            ['rapl.config', 'input', 'output'])
        self._remove_project_folder()

class TestControllerReadMapping(TestController):

    def test_read_mapping(self):
        self.controller.start_project()
        self.controller.paths._set_folder_names()
        self.controller.paths._set_static_file_names()
        self._generate_input_fasta_files()
        # If number of reads is less than the number of threads
        # segemehl stops. So set the number of threads to 1
        self.controller.parameters.segemehl_number_of_threads = 1
        self.controller.map_reads()
        self._remove_project_folder()

class TestControllerGRCreation(TestController):

    def test_create_gr_files(self):
        self.controller.start_project()
        self.controller.paths._set_folder_names()
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

class TestControllerAnnotationOverlap(TestController):
    
    def test_search_annotation_overlaps(self):
        self.controller.start_project()
        self.controller.paths._set_folder_names()
        self.controller.paths._set_static_file_names()
        self._generate_input_fasta_files()
        self._generate_mapping_files()
        self._generate_annotation_files()
        self.controller.search_annotation_overlaps()
        self.assertEqual(
            sorted(list(os.listdir(
                        self.controller.paths.annotation_hit_folder))),
            ['libbar.fa_annotation_overlaps.txt', 
             'libfoo.fa_annotation_overlaps.txt'])
        print(self.controller.paths.annotation_hit_folder)
        self.assertEqual(open(
                self.controller.paths.annotation_hit_folder + 
                "/" + 'libbar.fa_annotation_overlaps.txt').read(), 
                         self.example_data.overlap_output_1)
        self.assertEqual(open(
                self.controller.paths.annotation_hit_folder + 
                "/" + 'libfoo.fa_annotation_overlaps.txt').read(),
                         self.example_data.overlap_output_2)
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

if __name__ == "__main__":
    unittest.main()

