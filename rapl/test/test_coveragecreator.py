import sys
import os
import unittest
sys.path.append(".")
from libs.coveragecreator import CoverageCreator
from io import StringIO
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
        self.assertListEqual(
            sorted(self.coverage_creator.replicons_and_coverages[
                "reverse"].keys()),
            ['chrom', 'plasmid1', 'plasmid2'])
        self.assertListEqual(
            sorted(
                self.coverage_creator.replicons_and_coverages["forward"].keys()),
            ['chrom', 'plasmid1', 'plasmid2'])
        self.assertListEqual(
            self.coverage_creator.replicons_and_coverages["reverse"]["chrom"],
            [0] * 1500)
        self.assertListEqual(
            self.coverage_creator.replicons_and_coverages["forward"]["chrom"],
            [0] * 1500)
        self.assertListEqual(
            self.coverage_creator.replicons_and_coverages["reverse"]["plasmid1"],
            [0] * 100)
        self.assertListEqual(
            self.coverage_creator.replicons_and_coverages["forward"]["plasmid1"],
            [0] * 100)

    def test_count_coverage_1(self):
        """Check correct start at first list element"""
        self._generate_bam_file(
            self.example_data.sam_content_1, self._sam_bam_prefix)
        self.coverage_creator.init_coverage_lists("dummy.bam")
        self.coverage_creator.count_coverage("dummy.bam")
        self.assertEqual(
            len(self.coverage_creator.replicons_and_coverages[
                "forward"]["chrom"]),
            1500)
        self.assertListEqual(
            self.coverage_creator.replicons_and_coverages[
                "forward"]["chrom"][0:15],
            [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
             0.0, 0.0, 0.0, 0.0, 0.0])
        self.assertListEqual(
            self.coverage_creator.replicons_and_coverages[
                "reverse"]["chrom"][0:15],
            [-5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0,
             0.0, 0.0, 0.0, 0.0, 0.0])

    def test_count_coverage_2(self):
        """Consider how often a read is mapped. Mappings of reads that
        are aligned to several location contribute only fractions to
        the counting.

        """
        self._generate_bam_file(
            self.example_data.sam_content_2, self._sam_bam_prefix)
        self.coverage_creator.init_coverage_lists("dummy.bam")
        self.coverage_creator.count_coverage("dummy.bam")
        self.assertListEqual(
            self.coverage_creator.replicons_and_coverages[
                "forward"]["chrom"][0:15],
            [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
             0.0, 0.0, 0.0, 0.0, 0.0])

    def test_count_coverage_3(self):
        """If read_count_splitting is not set to True then every
        mapping is counted as one to each of the matching position
        independent how often its read is mapped in in total.
        """
        self._generate_bam_file(
            self.example_data.sam_content_2, self._sam_bam_prefix)
        self.coverage_creator.init_coverage_lists("dummy.bam")
        self.coverage_creator.count_coverage(
            "dummy.bam", read_count_splitting=False)
        self.assertListEqual(
            self.coverage_creator.replicons_and_coverages[
                "forward"]["chrom"][0:15],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
             0.0, 0.0, 0.0, 0.0, 0.0])

    def test_count_coverage_4(self):
        """If uniqueley_mapped_only is True scrip any mapping of read
        that are aligned to more than on location.
        """
        self._generate_bam_file(
            self.example_data.sam_content_3, self._sam_bam_prefix)
        self.coverage_creator.init_coverage_lists("dummy.bam")
        self.coverage_creator.count_coverage(
            "dummy.bam", uniqueley_mapped_only=True)
        self.assertListEqual(
            self.coverage_creator.replicons_and_coverages[
                "forward"]["chrom"][0:15],
            [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
             0.0, 0.0, 0.0, 0.0, 0.0])

    def test_write_to_wiggle_file_1(self):
        self._generate_bam_file(
            self.example_data.sam_content_3, self._sam_bam_prefix)
        self.coverage_creator.replicons_and_coverages = {
            "forward" : {"mychrom": [1.0, 2.0, 3.0]},
            "reverse" : {"mychrom": [1.0, 2.0, 3.0]}}
        output_fh = StringIO()
        self.coverage_creator._write_to_wiggle_file(
            output_fh, "boing", 1.0, "forward")
        self.assertEqual(
            output_fh.getvalue(),
            "track type=wiggle_0 name=\"boing_forward\"\nvariableStep "
            "chrom=mychrom span=1\n1 1.0\n2 2.0\n3 3.0\n")

    def test_write_to_wiggle_file_2(self):
        """Write wiggle file with factor."""
        self._generate_bam_file(
            self.example_data.sam_content_3, self._sam_bam_prefix)
        self.coverage_creator.replicons_and_coverages = {
            "forward" : {"mychrom": [1.0, 2.0, 3.0]},
            "reverse" : {"mychrom": [1.0, 2.0, 3.0]}}
        output_fh = StringIO()
        self.coverage_creator._write_to_wiggle_file(
            output_fh, "boing", 5.0, "forward")
        self.assertEqual(
            output_fh.getvalue(),
            "track type=wiggle_0 name=\"boing_forward\"\nvariableStep "
            "chrom=mychrom span=1\n1 5.0\n2 10.0\n3 15.0\n")

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
myread:005	6	chrom	5	255	10M	*	0	0	GTGGACAACC	*	NM:i:1	MD:Z:11T3	NH:i:9	XI:i:1	XA:Z:Q
"""

if __name__ == "__main__":
    unittest.main()
