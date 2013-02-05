import sys
import os
import unittest
sys.path.append(".")
from libs.coveragecreator import CoverageCreator
import pysam

class TestCoverageCreator(unittest.TestCase):

    def setUp(self):
        self.coverage_creator = CoverageCreator()
        self.example_data = ExampleData()
        self._sam_bam_prefix = "dummy"

    def tearDown(self):
        for suffix in [".sam", ".bam", ".bam.bai"]:
            os.remove(self._sam_bam_prefix + suffix)

    def _generate_bam_file(self, sam_content, file_prefix):
        sam_file = "%s.sam" % file_prefix
        bam_file = "%s.bam" % file_prefix
        sam_fh = open(sam_file, "w")
        sam_fh.write(sam_content)
        sam_fh.close()
        pysam.view("-Sb", "-o%s" % bam_file, sam_file)
        pysam.index(bam_file)

    def test_init_coverage_lists(self):
        self._generate_bam_file(
            self.example_data.sam_content_1, self._sam_bam_prefix)
        self.coverage_creator.init_coverage_lists("dummy.bam")
        self.assertEqual(
            sorted(self.coverage_creator.elements_and_coverages["minus"].keys()),
            ['chrom', 'plasmid1', 'plasmid2'])
        self.assertEqual(
            sorted(self.coverage_creator.elements_and_coverages["plus"].keys()),
            ['chrom', 'plasmid1', 'plasmid2'])
        self.assertEqual(
            self.coverage_creator.elements_and_coverages["minus"]["chrom"],
            [0] * 1500)
        self.assertEqual(
            self.coverage_creator.elements_and_coverages["plus"]["chrom"],
            [0] * 1500)
        self.assertEqual(
            self.coverage_creator.elements_and_coverages["minus"]["plasmid1"],
            [0] * 100)
        self.assertEqual(
            self.coverage_creator.elements_and_coverages["plus"]["plasmid1"],
            [0] * 100)

    def test_count_coverage_1(self):
        self._generate_bam_file(
            self.example_data.sam_content_1, self._sam_bam_prefix)
        self.coverage_creator.init_coverage_lists("dummy.bam")
        self.coverage_creator.count_coverage("dummy.bam")
        self.assertEqual(
            len(self.coverage_creator.elements_and_coverages["plus"]["chrom"]),
            1500)
        # Check correct start at first list element
        self.assertEqual(
            self.coverage_creator.elements_and_coverages["plus"]["chrom"][0:15],
            [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
             0.0, 0.0, 0.0, 0.0, 0.0])

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
"""

    sam_content_2 = """@HD	VN:1.0
@SQ	SN:chrom	LN:1500
@SQ	SN:plasmid1	LN:100
@SQ	SN:plasmid2	LN:200
myread:2306	0	chrom	1	255	13M	*	0	0	TAACTGGCTGTGG	*	NM:i:0	MD:Z:13	NH:i:1	XI:i:0	XA:Z:Q
myread:2107	0	chrom	10	255	15M	*	0	0	GTGGACAACCGATTT	*	NM:i:1	MD:Z:11T3	NH:i:1	XI:i:1	XA:Z:Q
myread:2304	0	chrom	20	255	23M	*	0	0	GCTGTCGAATAGAAATATAGCTG	*	NM:i:1	MD:Z:0C22	NH:i:1	XI:i:0	XA:Z:Q
myread:1105	0	chrom	50	255	19M	*	0	0	CCTGTCGAATAGAAATATA	*	NM:i:0	MD:Z:19	NH:i:1	XI:i:0	XA:Z:Q
myread:1207	0	chrom	60	255	24M	*	0	0	CCTGTCGAATAGAAATATAGCTGG	*	NM:i:0	MD:Z:24	NH:i:1	XI:i:0	XA:Z:Q
myread:2305	0	chrom	80	255	24M	*	0	0	ACTGTCGAATAGAAATATAGCTGG	*	NM:i:1	MD:Z:0C23	NH:i:1	XI:i:0	XA:Z:Q
"""

if __name__ == "__main__":
    unittest.main()

