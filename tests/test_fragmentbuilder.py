import sys
import os
import unittest

sys.path.append(".")
from reademptionlib.fragmentbuilder import FragmentBuilder
import pysam


class TestFragmentBuilder(unittest.TestCase):
    def setUp(self):
        #self.species_references = ["chrom", "plasmid1", "plasmid2"]
        self.max_fragment_length = False
        self.fragment_builder = FragmentBuilder(self.max_fragment_length)
        self.example_data = ExampleData()

        self._no_fragments = "no_fragments"
        self.input_bam_no_fragments_path = "no_fragments.bam"

        self._with_fragments_expected = "with_fragments_expected"
        self.bam_with_fragments_expected_path = "with_fragments_expected.bam"
        self.sam_with_fragments_expected_path = "with_fragments_expected.sam"

        #self._with_fragments_expected = "with_fragments_expected"
        self._with_fragments_built = "with_fragments_built"
        self.bam_with_fragments_built_path = "with_fragments_built.bam"
        self.sam_with_fragments_built_path = "with_fragments_built.sam"

    def tearDown(self):
        for suffix in [".sam", ".bam", ".bam.bai"]:
            for path in [self._no_fragments, self._with_fragments_expected, self._with_fragments_built]:
                if os.path.exists(path + suffix) is True:
                    os.remove(path + suffix)

    def _generate_bam_and_sam_file(self, sam_content, file_prefix):
        sam_file = "{}.sam".format(file_prefix)
        bam_file = "{}.bam".format(file_prefix)
        sam_fh = open(sam_file, "w")
        sam_fh.write(sam_content)
        sam_fh.close()
        pysam.view("-h", "-o", bam_file, sam_file, catch_stdout=False)
        pysam.index(bam_file)


    def test_fragment_builder(self):
        # generate input bam and sam file without fragments
        self._generate_bam_and_sam_file(
            self.example_data.sam_no_fragments, self._no_fragments
        )
        # generate expected output bam and sam file with fragments
        self._generate_bam_and_sam_file(
            self.example_data.sam_with_fragments_expected, self._with_fragments_expected
        )
        # Run the fragment builder with the input bam
        self.fragment_builder.build_bam_file_with_fragments(self.input_bam_no_fragments_path, self.bam_with_fragments_built_path)

        # Create a sam file of the bam file that was created by the fragment builder
        pysam.view("-h", "-o", self.sam_with_fragments_built_path, self.bam_with_fragments_built_path, catch_stdout=False)

        # Compare the expected sam file content with the created sam file content
        # Can be used to compare the content of two files line wise:
        with open(self.sam_with_fragments_expected_path, 'r') as expected_sam, open(self.sam_with_fragments_built_path, 'r') as calculated_sam:
            expected_contend = expected_sam.readlines()
            calculated_contend = calculated_sam.readlines()
            for l1, l2 in zip(expected_contend, calculated_contend):
                if (l1.startswith("@") or l2.startswith("@")):
                    continue
                if l1 != l2:
                    print(f"expected line: '{l1}' does not match calculated line '{l2}'")
            assert expected_contend == calculated_contend

class ExampleData(object):
    sam_no_fragments = """@HD	VN:1.0	SO:coordinate
@SQ	SN:chr1	LN:1140
@SQ	SN:NC_007795.1	LN:3360
@RG	ID:A1	SM:sample1	LB:library1	PU:unit1	PL:illumina
@PG	ID:segemehl	VN:0.3.4	CL:segemehl.x --query reademption_analysis_dual_paired_end/output/align/processed_reads/library_one_p1_processed.fa.gz --mate reademption_analysis_dual_paired_end/output/align/processed_reads/library_one_p2_processed.fa.gz --index reademption_analysis_dual_paired_end/output/align/index/index.idx --database reademption_analysis_dual_paired_end/input/human_reference_sequences/GRCh38.p10.genome_short_for_paired_end.fa reademption_analysis_dual_paired_end/input/staphylococcus_reference_sequences/GCF_000013425.1_ASM1342v1_genomic_short_for_paired_end.fa --outfile reademption_analysis_dual_paired_end/output/align/alignments/library_one_alignments_final.bam --bamabafixoida --hitstrategy 1 --accuracy 100 --evalue 5.0 --threads 1 --splits --nomatchfilename reademption_analysis_dual_paired_end/output/align/unaligned_reads/library_one_unaligned.fa
@PG	ID:samtools	PN:samtools	PP:segemehl	VN:1.10 (pysam)	CL:samtools view -b -F 4 -o reademption_analysis_dual_paired_end/output/align/alignments/library_one_alignments_final.bam_filtered reademption_analysis_dual_paired_end/output/align/alignments/library_one_alignments_final.bam
@PG	ID:samtools.1	PN:samtools	PP:samtools	VN:1.10 (pysam)	CL:samtools sort -o reademption_analysis_dual_paired_end/output/align/alignments/library_one_alignments_final.bam_sorted reademption_analysis_dual_paired_end/output/align/alignments/library_one_alignments_final.bam
SN7001299:308:CAPEHACXX:5:2301:2152:1	355	chr1	1	1	20=	=	41	60	ACCCTAACCCTAACCCTAAC	*	HI:i:1	NH:i:2	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:1	147	chr1	41	1	20=	=	541	520	TAACCCTAACCCTAACCCTA	*	HI:i:0	NH:i:2	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:1	403	chr1	41	1	20=	=	1	-60	TAACCCTAACCCTAACCCTA	*	HI:i:1	NH:i:2	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:2	99	chr1	61	1	20=	=	101	60	TCTGACCTGAGGAGAACTGT	*	HI:i:0	NH:i:4	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:2	355	chr1	61	1	20=	=	221	180	TCTGACCTGAGGAGAACTGT	*	HI:i:1	NH:i:4	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:2	147	chr1	101	1	20=	=	61	-60	CCGAAATCTGTGCAGAGGAC	*	HI:i:0	NH:i:4	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:2	403	chr1	101	1	20=	=	181	100	CCGAAATCTGTGCAGAGGAC	*	HI:i:2	NH:i:4	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:2	355	chr1	181	1	20=	=	101	-100	TCTGACCTGAGGAGAACTGT	*	HI:i:2	NH:i:4	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:2	355	chr1	181	1	20=	=	221	60	TCTGACCTGAGGAGAACTGT	*	HI:i:3	NH:i:4	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:2	403	chr1	221	1	20=	=	61	-180	CCGAAATCTGTGCAGAGGAC	*	HI:i:1	NH:i:4	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:2	403	chr1	221	1	20=	=	181	-60	CCGAAATCTGTGCAGAGGAC	*	HI:i:3	NH:i:4	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:3	163	chr1	301	1	20=	=	341	60	GCGGTCATGCGCCCGTATAT	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:3	83	chr1	341	1	20=	=	301	-60	TTTATTATGCGGCATTAGCG	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:4	89	chr1	401	1	20=	*	0	0	CGCGCCGGTTCAGGTGCAGA	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:5	73	chr1	481	1	20=	*	0	0	CGCAGGCGCAGAGACACATG	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:1	99	chr1	541	1	20=	=	41	-520	ACCCTAACCCTAACCCTAAC	*	HI:i:0	NH:i:2	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:6	339	chr1	641	1	20=	=	661	40	GTCAGCTGCAATCTATGCCG	*	HI:i:1	NH:i:2	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:6	163	chr1	661	1	20=	=	701	60	TGACTCCGCGCGTACACGTC	*	HI:i:0	NH:i:2	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:6	419	chr1	661	1	20=	=	641	-40	TGACTCCGCGCGTACACGTC	*	HI:i:1	NH:i:2	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:6	83	chr1	701	1	20=	=	661	-60	GTCAGCTGCAATCTATGCCG	*	HI:i:0	NH:i:2	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:7	99	chr1	731	1	20=	=	741	30	CATTGAGATCTGACGGAGAG	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:7	147	chr1	741	1	20=	=	731	-30	TGACGGAGAGCACCTCGAGC	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:8	99	chr1	781	1	20=	=	781	-20	GCTATAATATGCAATGTGCG	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:8	147	chr1	781	1	20=	=	781	20	GCTATAATATGCAATGTGCG	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:9	147	chr1	846	1	20=	=	851	25	CGTCACGTATTCAGTCGATG	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:9	99	chr1	851	1	20=	=	846	-25	CGTATTCAGTCGATGTGAGA	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:10	99	chr1	901	1	43=7S	=	1001	120	AAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAGCTTAAAA	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:43	XX:i:1	XY:i:50	XI:i:0	XH:i:1	XJ:i:1	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:10	147	chr1	1001	1	20=	=	901	-120	TAGCCTAGCTTAGCGTAGCT	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:11	163	chr1	1021	1	20=	=	1031	30	CGATCGGCGTATTAGTCGCG	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:11	83	chr1	1031	1	20=	=	1021	-30	ATTAGTCGCGAGTTATCGAG	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:12	83	chr1	1086	1	20=	=	1091	25	ATGCACGTGTGTGACACAGT	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:12	163	chr1	1091	1	20=	=	1086	-25	CGTGTGTGACACAGTGACGG	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
"""

    sam_with_fragments_expected="""@HD	VN:1.0	SO:coordinate
@SQ	SN:chr1	LN:1140
@SQ	SN:NC_007795.1	LN:3360
@RG	ID:A1	SM:sample1	LB:library1	PU:unit1	PL:illumina
@PG	ID:segemehl	VN:0.3.4	CL:segemehl.x --query reademption_analysis_dual_paired_end/output/align/processed_reads/library_one_p1_processed.fa.gz --mate reademption_analysis_dual_paired_end/output/align/processed_reads/library_one_p2_processed.fa.gz --index reademption_analysis_dual_paired_end/output/align/index/index.idx --database reademption_analysis_dual_paired_end/input/human_reference_sequences/GRCh38.p10.genome_short_for_paired_end.fa reademption_analysis_dual_paired_end/input/staphylococcus_reference_sequences/GCF_000013425.1_ASM1342v1_genomic_short_for_paired_end.fa --outfile reademption_analysis_dual_paired_end/output/align/alignments/library_one_alignments_final.bam --bamabafixoida --hitstrategy 1 --accuracy 100 --evalue 5.0 --threads 1 --splits --nomatchfilename reademption_analysis_dual_paired_end/output/align/unaligned_reads/library_one_unaligned.fa
@PG	ID:samtools	PN:samtools	PP:segemehl	VN:1.10 (pysam)	CL:samtools view -b -F 4 -o reademption_analysis_dual_paired_end/output/align/alignments/library_one_alignments_final.bam_filtered reademption_analysis_dual_paired_end/output/align/alignments/library_one_alignments_final.bam
@PG	ID:samtools.1	PN:samtools	PP:samtools	VN:1.10 (pysam)	CL:samtools sort -o reademption_analysis_dual_paired_end/output/align/alignments/library_one_alignments_final.bam_sorted reademption_analysis_dual_paired_end/output/align/alignments/library_one_alignments_final.bam
@PG	ID:samtools.2	PN:samtools	PP:samtools.1	VN:1.10 (pysam)	CL:samtools view -h -o no_fragments.bam no_fragments.sam
@PG	ID:samtools.3	PN:samtools	PP:samtools.2	VN:1.10 (pysam)	CL:samtools sort -o with_fragments_built.bam_sorted with_fragments_built.bam
@PG	ID:samtools.4	PN:samtools	PP:samtools.3	VN:1.10 (pysam)	CL:samtools view -h -o with_fragments_built.sam with_fragments_built.bam
SN7001299:308:CAPEHACXX:5:2301:2152:1	355	chr1	1	1	60=	=	41	60	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	*	HI:i:1	NH:i:2	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:1	99	chr1	41	1	520=	=	41	-520	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	*	HI:i:0	NH:i:2	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:2	99	chr1	61	1	60=	=	101	60	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	*	HI:i:0	NH:i:4	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:2	355	chr1	61	1	180=	=	221	180	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	*	HI:i:1	NH:i:4	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:2	355	chr1	101	1	100=	=	101	-100	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	*	HI:i:2	NH:i:4	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:2	355	chr1	181	1	60=	=	221	60	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	*	HI:i:3	NH:i:4	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:3	83	chr1	301	1	60=	=	301	-60	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:4	89	chr1	401	1	20=	*	0	0	CGCGCCGGTTCAGGTGCAGA	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:5	73	chr1	481	1	20=	*	0	0	CGCAGGCGCAGAGACACATG	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:6	339	chr1	641	1	40=	=	661	40	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	*	HI:i:1	NH:i:2	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:6	83	chr1	661	1	60=	=	661	-60	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	*	HI:i:0	NH:i:2	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:7	99	chr1	731	1	30=	=	741	30	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:8	99	chr1	781	1	20=	=	781	-20	NNNNNNNNNNNNNNNNNNNN	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:9	99	chr1	851	1	15=	=	846	-25	NNNNNNNNNNNNNNN	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:10	99	chr1	901	1	120=	=	1001	120	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:43	XX:i:1	XY:i:50	XI:i:0	XH:i:1	XJ:i:1	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:11	83	chr1	1021	1	30=	=	1021	-30	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
SN7001299:308:CAPEHACXX:5:2301:2152:12	83	chr1	1091	1	15=	=	1091	25	NNNNNNNNNNNNNNN	*	HI:i:0	NH:i:1	NM:i:0	MD:Z:20	RG:Z:A1	YZ:Z:0
"""


if __name__ == "__main__":
    unittest.main()
