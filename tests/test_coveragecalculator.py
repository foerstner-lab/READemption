import sys
import os
import unittest

sys.path.append(".")
from reademptionlib.coveragecalculator import CoverageCalculator
from io import StringIO
import pysam


class TestCoverageCalculator(unittest.TestCase):
    def setUp(self):
        self.species_references = ["chrom", "plasmid1", "plasmid2"]
        self.coverage_calculator = CoverageCalculator(self.species_references, count_cross_aligned_reads=True)
        self.example_data = ExampleData()
        self._sam_bam_prefix = "dummy"

    def tearDown(self):
        for suffix in [".sam", ".bam", ".bam.bai"]:
            if os.path.exists(self._sam_bam_prefix + suffix) is True:
                os.remove(self._sam_bam_prefix + suffix)

    def _generate_bam_file(self, sam_content, file_prefix):
        sam_file = "{}.sam".format(file_prefix)
        bam_file = "{}.bam".format(file_prefix)
        sam_fh = open(sam_file, "w")
        sam_fh.write(sam_content)
        sam_fh.close()
        pysam.view("-Sb", "-o{}".format(bam_file), sam_file, catch_stdout=False)
        pysam.index(bam_file)
        self._bam = pysam.Samfile(bam_file)
        #self._bam_path = bam_file

    def test_init_coverage_list(self):
        self.coverage_calculator._init_coverage_list(10)
        self.assertListEqual(
            sorted(self.coverage_calculator._coverages.keys()),
            ["forward", "reverse"],
        )
        self.assertListEqual(
            self.coverage_calculator._coverages["forward"].tolist(), [0.0] * 10
        )
        self.assertListEqual(
            self.coverage_calculator._coverages["reverse"].tolist(), [0.0] * 10
        )

    def test_calc_coverage_1(self):
        """Check correct start at first list element"""
        self._generate_bam_file(
            self.example_data.sam_content_1, self._sam_bam_prefix
        )
        self.coverage_calculator._init_coverage_list(self._bam.lengths[0])
        self.coverage_calculator._calc_coverage("chrom", self._bam)
        self.assertEqual(
            len(self.coverage_calculator._coverages["forward"]), 1500
        )
        self.assertListEqual(
            self.coverage_calculator._coverages["forward"][0:15].tolist(),
            [
                5.0,
                5.0,
                5.0,
                5.0,
                5.0,
                5.0,
                5.0,
                5.0,
                5.0,
                5.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
        )
        self.assertListEqual(
            self.coverage_calculator._coverages["reverse"][0:15].tolist(),
            [
                -5.0,
                -5.0,
                -5.0,
                -5.0,
                -5.0,
                -5.0,
                -5.0,
                -5.0,
                -5.0,
                -5.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
        )

    def test_calc_coverage_2(self):
        """Consider how often a read is mapped. Mappings of reads that
        are aligned to several location contribute only fractions to
        the counting.

        """
        self._generate_bam_file(
            self.example_data.sam_content_2, self._sam_bam_prefix
        )
        self.coverage_calculator._init_coverage_list(self._bam.lengths[0])
        self.coverage_calculator._calc_coverage("chrom", self._bam)
        self.assertListEqual(
            self.coverage_calculator._coverages["forward"][0:15].tolist(),
            [
                0.5,
                0.5,
                0.5,
                0.5,
                0.5,
                0.5,
                0.5,
                0.5,
                0.5,
                0.5,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
        )

    def test_calc_coverage_3(self):
        """If read_count_splitting is set to False then every
        mapping is counted as one to each of the matching position
        independent how often its read is mapped in in total.
        """
        self.coverage_calculator = CoverageCalculator(self.species_references,
            read_count_splitting=False, count_cross_aligned_reads=True
        )
        self._generate_bam_file(
            self.example_data.sam_content_2, self._sam_bam_prefix
        )
        self.coverage_calculator._init_coverage_list(self._bam.lengths[0])
        self.coverage_calculator._calc_coverage("chrom", self._bam)
        self.assertListEqual(
            self.coverage_calculator._coverages["forward"][0:15].tolist(),
            [
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
        )

    def test_calc_coverage_4(self):
        """If uniqueley_aligned_only is True skip any mapping of read
        that are aligned to more than on location.
        """
        self.coverage_calculator = CoverageCalculator(self.species_references,
            uniquely_aligned_only=True, count_cross_aligned_reads=True
        )
        self._generate_bam_file(
            self.example_data.sam_content_3, self._sam_bam_prefix
        )
        self.coverage_calculator._init_coverage_list(self._bam.lengths[0])
        self.coverage_calculator._calc_coverage("chrom", self._bam)
        self.assertListEqual(
            self.coverage_calculator._coverages["forward"][0:15].tolist(),
            [
                3.0,
                3.0,
                3.0,
                3.0,
                3.0,
                3.0,
                3.0,
                3.0,
                3.0,
                3.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
        )

    def test_calc_coverage_5(self):
        """If first_base_only is True only the first nucleotide of a
        mapping is considered.
        """
        self.coverage_calculator = CoverageCalculator(self.species_references,
            coverage_style="first_base_only", count_cross_aligned_reads=True
        )
        self._generate_bam_file(
            self.example_data.sam_content_1, self._sam_bam_prefix
        )
        self.coverage_calculator._init_coverage_list(self._bam.lengths[0])
        self.coverage_calculator._calc_coverage("chrom", self._bam)
        self.assertListEqual(
            self.coverage_calculator._coverages["forward"][0:15].tolist(),
            [
                5.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
        )
        self.assertListEqual(
            self.coverage_calculator._coverages["reverse"][0:15].tolist(),
            [
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                -5.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
        )


class ExampleData(object):

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


if __name__ == "__main__":
    unittest.main()
