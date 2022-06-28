import unittest
import os
import sys

sys.path.append(".")
from reademptionlib.genewisequanti import GeneWiseQuantification
import pysam


class Gff3EntryMoc(object):
    def __init__(self, seq_id, start, end, strand):
        self.seq_id = seq_id
        self.start = start
        self.end = end
        self.strand = strand


class TestGeneWiseQuantification(unittest.TestCase):
    def setUp(self):
        self.example_data = ExampleData()
        self._sam_bam_prefix = "dummy"
        self.gene_wise_quantification = GeneWiseQuantification("references_by_species_dummy")

    def tearDown(self):
        for suffix in [".sam", ".bam", ".bam.bai"]:
            os.remove(self._sam_bam_prefix + suffix)

    def test_overlapping_alignments_1(self):
        self._generate_bam_file(
            self.example_data.sam_content_1, self._sam_bam_prefix
        )
        self.gene_wise_quantification = GeneWiseQuantification("references_by_species_dummy")
        sam = pysam.Samfile(self._sam_bam_prefix + ".bam")
        # Overlap with all mappings on the forward strand
        self.assertListEqual(
            self._mapping_ids(
                self.gene_wise_quantification._overlapping_alignments(
                    sam, Gff3EntryMoc("chrom", 1, 100, "+")
                )
            ),
            ["myread:01", "myread:02", "myread:03", "myread:04", "myread:05"],
        )
        # Overlapping with no mapping
        self.assertListEqual(
            self._mapping_ids(
                self.gene_wise_quantification._overlapping_alignments(
                    sam, Gff3EntryMoc("chrom", 1, 5, "+")
                )
            ),
            [],
        )
        # Overlapping by 1 based - in the 5' end of the reads
        self.assertListEqual(
            self._mapping_ids(
                self.gene_wise_quantification._overlapping_alignments(
                    sam, Gff3EntryMoc("chrom", 1, 10, "+")
                )
            ),
            ["myread:01", "myread:02", "myread:03", "myread:04", "myread:05"],
        )
        # No overlap - gene very close upstream of the reads
        self.assertListEqual(
            self._mapping_ids(
                self.gene_wise_quantification._overlapping_alignments(
                    sam, Gff3EntryMoc("chrom", 1, 9, "+")
                )
            ),
            [],
        )
        # Overlapping by 1 based - in the 3' end of the reads
        self.assertListEqual(
            self._mapping_ids(
                self.gene_wise_quantification._overlapping_alignments(
                    sam, Gff3EntryMoc("chrom", 19, 23, "+")
                )
            ),
            ["myread:01", "myread:02", "myread:03", "myread:04", "myread:05"],
        )
        # Overlapping by 1 based - in the 3' end of the reads but on the wrong strand
        self.assertListEqual(
            self._mapping_ids(
                self.gene_wise_quantification._overlapping_alignments(
                    sam, Gff3EntryMoc("chrom", 19, 23, "-")
                )
            ),
            [],
        )
        # No overlap - very close downstream of the reads
        self.assertListEqual(
            self._mapping_ids(
                self.gene_wise_quantification._overlapping_alignments(
                    sam, Gff3EntryMoc("chrom", 20, 23, "+")
                )
            ),
            [],
        )
        # Overlapping by 1 based - in the 3' end of the reads on the opposite strand, without strand specificity
        self.gene_wise_quantification._strand_specific = False
        self.assertListEqual(
            self._mapping_ids(
                self.gene_wise_quantification._overlapping_alignments(
                    sam, Gff3EntryMoc("chrom", 19, 23, "-")
                )
            ),
            ["myread:01", "myread:02", "myread:03", "myread:04", "myread:05"],
        )
        # Overlapping by 1 based - in the 3' end of the reads on the opposite strand, with strand specificity
        # only allowing antisense overlaps
        self.gene_wise_quantification._strand_specific = True
        self.gene_wise_quantification._antisense_only = True
        self.assertListEqual(
            self._mapping_ids(
                self.gene_wise_quantification._overlapping_alignments(
                    sam, Gff3EntryMoc("chrom", 19, 100, "-")
                )
            ),
            ["myread:01", "myread:02", "myread:03", "myread:04", "myread:05"],
        )

    def test_overlapping_alignments_2(self):
        """Extraction of overlapping reads - with a non-default
        minimal overlap.
        """
        self._generate_bam_file(
            self.example_data.sam_content_1, self._sam_bam_prefix
        )
        self.gene_wise_quantification = GeneWiseQuantification("references_by_species_dummy")
        self.gene_wise_quantification._min_overlap = 5
        sam = pysam.Samfile(self._sam_bam_prefix + ".bam")
        # 1 overlapping base in the 5' end of the reads => not enough
        self.assertListEqual(
            self._mapping_ids(
                self.gene_wise_quantification._overlapping_alignments(
                    sam, Gff3EntryMoc("chrom", 1, 10, "+")
                )
            ),
            [],
        )
        # 4 overlapping base in the 5' end of the reads => not enough
        self.assertListEqual(
            self._mapping_ids(
                self.gene_wise_quantification._overlapping_alignments(
                    sam, Gff3EntryMoc("chrom", 1, 13, "+")
                )
            ),
            [],
        )
        # 5 overlapping base in the 5' end of the reads => okay
        self.assertListEqual(
            self._mapping_ids(
                self.gene_wise_quantification._overlapping_alignments(
                    sam, Gff3EntryMoc("chrom", 1, 14, "+")
                )
            ),
            ["myread:01", "myread:02", "myread:03", "myread:04", "myread:05"],
        )
        # 1 overlapping base in the 3' end of the reads => not enough
        self.assertListEqual(
            self._mapping_ids(
                self.gene_wise_quantification._overlapping_alignments(
                    sam, Gff3EntryMoc("chrom", 19, 23, "+")
                )
            ),
            [],
        )
        # 4 overlapping base in the 3' end of the reads => not enough
        self.assertListEqual(
            self._mapping_ids(
                self.gene_wise_quantification._overlapping_alignments(
                    sam, Gff3EntryMoc("chrom", 16, 23, "+")
                )
            ),
            [],
        )
        # 5 overlapping base in the 3' end of the reads => okay
        self.assertListEqual(
            self._mapping_ids(
                self.gene_wise_quantification._overlapping_alignments(
                    sam, Gff3EntryMoc("chrom", 15, 23, "+")
                )
            ),
            ["myread:01", "myread:02", "myread:03", "myread:04", "myread:05"],
        )

    def _mapping_ids(self, mappings):
        return [mapping.qname for mapping in mappings]

    def _generate_bam_file(self, sam_content, file_prefix):
        sam_file = "{}.sam".format(file_prefix)
        bam_file = "{}.bam".format(file_prefix)
        sam_fh = open(sam_file, "w")
        sam_fh.write(sam_content)
        sam_fh.close()
        pysam.view("-Sb", "-o{}".format(bam_file), sam_file, catch_stdout=False)
        pysam.index(bam_file)


class ExampleData(object):

    sam_content_1 = """@HD	VN:1.0
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


if __name__ == "__main__":
    unittest.main()
