from reademptionlib.genewisequanti import GeneWiseQuantification


def data_gene_wise_quanti():
    gene_wise_quantification = GeneWiseQuantification()
    sam_bam_prefix = "dummy"
    sam_content = """@HD	VN:1.0
@SQ	SN:chrom	LN:1500
@SQ	SN:plasmid1	LN:100
@SQ	SN:plasmid2	LN:200
myread:01	0	chrom	10	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:02	0	chrom	10	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:03	0	chrom	10	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:04	0	chrom	10	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:05	0	chrom	10	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:06	16	chrom	35	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:07	16	chrom	35	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:08	16	chrom	35	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:09	16	chrom	35	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:10	16	chrom	35	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
"""

    global gene_wise_quantification
    global sam_bam_prefix
    global sam_content
    
