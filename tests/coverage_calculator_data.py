from reademptionlib.coveragecalculator import CoverageCalculator


def data_coverage_calculator():
    coverage_calculator = CoverageCalculator()
    sam_bam_prefix = "dummy"
    sam_content_1 = """@HD	VN:1.0
@SQ	SN:chrom	LN:1500
@SQ	SN:plasmid1	LN:100
@SQ	SN:plasmid2	LN:200
myread:001	0	chrom	1	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:002	0	chrom	1	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:003	0	chrom	1	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:004	0	chrom	1	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:005	0	chrom	1	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:006	16	chrom	1	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:007	16	chrom	1	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:008	16	chrom	1	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:009	16	chrom	1	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:010	16	chrom	1	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
"""
    sam_content_2 = """@HD	VN:1.0
@SQ	SN:chrom	LN:1500
@SQ	SN:plasmid1	LN:100
@SQ	SN:plasmid2	LN:200
myread:001	0	chrom	1	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:2	XI:i:1	XA:Z:Q
"""

    sam_content_3 = """@HD	VN:1.0
@SQ	SN:chrom	LN:1500
@SQ	SN:plasmid1	LN:100
@SQ	SN:plasmid2	LN:200
myread:001	0	chrom	1	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:002	0	chrom	1	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:003	0	chrom	1	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:004	0	chrom	5	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:9	XI:i:1	XA:Z:Q
myread:005	0	chrom	5	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:9	XI:i:1	XA:Z:Q
myread:006	0	chrom	5	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:9	XI:i:1	XA:Z:Q
"""
    global coverage_calculator
    global sam_bam_prefix
    global sam_content_1
    global sam_content_2
    global sam_content_3
