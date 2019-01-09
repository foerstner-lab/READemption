import sys
import os
sys.path.append("./tests")
import coverage_calculator_data as ccd
from reademptionlib.coveragecalculator import CoverageCalculator
import pysam
import numpy as np


def setup_function(function):
    ccd.data_coverage_calculator()

    
def teardown_function(function):
    for suffix in [".sam", ".bam", ".bam.bai"]:
            if os.path.exists(ccd.sam_bam_prefix + suffix) is True:
                os.remove(ccd.sam_bam_prefix + suffix)


def test_init_coverage_list():
    ccd.coverage_calculator._init_coverage_list(10)
    assert sorted(ccd.coverage_calculator._coverages.keys()) == [
        "forward", "reverse"]
    assert (ccd.coverage_calculator._coverages["forward"]).all() == (
        np.array([0.0] * 10)).all()
    assert (ccd.coverage_calculator._coverages["reverse"]).all() == (
        np.array([0.0] * 10)).all()

    
def test_calc_coverage_1():
    """Check correct start at first list element"""
    bam_file_1 = generate_bam_file(ccd.sam_content_1, ccd.sam_bam_prefix)
    ccd.coverage_calculator._init_coverage_list(bam_file_1.lengths[0])
    ccd.coverage_calculator._calc_coverage("chrom", bam_file_1)
    assert len(ccd.coverage_calculator._coverages["forward"]) == 1500
    assert (ccd.coverage_calculator._coverages["reverse"][0:15]).all() == (
        np.array([-5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0,
                  -5.0, 0.0, 0.0, 0.0, 0.0, 0.0])).all()


def test_calc_coverage_2():
    """Consider how often a read is mapped. Mappings of reads that
        are aligned to several location contribute only fractions to
        the counting.
    """
    bam_file_2 = generate_bam_file(ccd.sam_content_2, ccd.sam_bam_prefix)
    ccd.coverage_calculator._init_coverage_list(bam_file_2.lengths[0])
    ccd.coverage_calculator._calc_coverage("chrom", bam_file_2)
    assert (ccd.coverage_calculator._coverages["forward"][0:15]).all() == (
        np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                  0.0, 0.0, 0.0, 0.0, 0.0])).all()


def test_calc_coverage_3():
    """If read_count_splitting is set to False then every
        mapping is counted as one to each of the matching position
        independent how often its read is mapped in in total.
    """
    coverage_calculator = CoverageCalculator(read_count_splitting=False)
    bam_file_3 = generate_bam_file(ccd.sam_content_2, ccd.sam_bam_prefix)
    coverage_calculator._init_coverage_list(bam_file_3.lengths[0])
    coverage_calculator._calc_coverage("chrom", bam_file_3)
    assert (coverage_calculator._coverages["forward"][0:15]).all() == (
        np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                  0.0, 0.0, 0.0, 0.0, 0.0])).all()


def test_calc_coverage_4():
    """If uniqueley_aligned_only is True skip any mapping of read
        that are aligned to more than on location.
    """
    coverage_calculator = CoverageCalculator(uniquely_aligned_only=True)
    bam_file_4 = generate_bam_file(ccd.sam_content_3, ccd.sam_bam_prefix)
    coverage_calculator._init_coverage_list(bam_file_4.lengths[0])
    coverage_calculator._calc_coverage("chrom", bam_file_4)
    assert (coverage_calculator._coverages["forward"][0:15]).all() == (
        np.array([3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
                  0.0, 0.0, 0.0, 0.0, 0.0])).all()


def test_calc_coverage_5():
    """If first_base_only is True only the first nucleotide of a
        mapping is considered.
    """
    coverage_calculator = CoverageCalculator(coverage_style="first_base_only")
    bam_file_5 = generate_bam_file(ccd.sam_content_1, ccd.sam_bam_prefix)
    coverage_calculator._init_coverage_list(bam_file_5.lengths[0])
    coverage_calculator._calc_coverage("chrom", bam_file_5)
    assert (coverage_calculator._coverages["forward"][0:15]).all() == (
        np.array([5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0])).all()
    assert (coverage_calculator._coverages["reverse"][0:15]).all() == (
        np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -5.0,
                  0.0, 0.0, 0.0, 0.0, 0.0])).all()


def generate_bam_file(sam_content, file_prefix):
    sam_file = "{}.sam".format(file_prefix)
    bam_file = "{}.bam".format(file_prefix)
    sam_fh = open(sam_file, "w")
    sam_fh.write(sam_content)
    sam_fh.close()
    pysam.view("-Sb", "-o{}".format(bam_file), sam_file, catch_stdout=False)
    pysam.index(bam_file)
    bam = pysam.Samfile(bam_file)
    return bam
