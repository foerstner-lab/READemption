import os
import sys
sys.path.append("./tests")
import gene_wise_quanti_data as gqd
import pysam


def setup_function(function):
    gqd.data_gene_wise_quanti()


def teardown_function(function):
    for suffix in [".sam", ".bam", ".bam.bai"]:
        os.remove(gqd.sam_bam_prefix + suffix)


def test_overlapping_alignments():
    generate_bam_file(gqd.sam_content, gqd.sam_bam_prefix)
    sam = pysam.Samfile(gqd.sam_bam_prefix + ".bam")
    # Overlap with all mappings
    assert mapping_ids(gqd.gene_wise_quantification._overlapping_alignments(
        sam, Gff3EntryMoc("chrom", 1, 100))) == [
            "myread:01", "myread:02", "myread:03", "myread:04", "myread:05",
            "myread:06", "myread:07", "myread:08", "myread:09", "myread:10"]
    # Overlapping with no mapping
    assert mapping_ids(gqd.gene_wise_quantification._overlapping_alignments(
        sam, Gff3EntryMoc("chrom", 1, 5))) == []
    # Overlapping by 1 based - in the 5' end of the reads
    assert mapping_ids(gqd.gene_wise_quantification._overlapping_alignments(
        sam, Gff3EntryMoc("chrom", 1, 10))) == [
            "myread:01", "myread:02", "myread:03", "myread:04", "myread:05"]
    # No overlap - gene very close upstream of the reads
    assert mapping_ids(gqd.gene_wise_quantification._overlapping_alignments(
        sam, Gff3EntryMoc("chrom", 1, 9))) == []
    # Overlapping by 1 based - in the 3' end of the reads
    assert mapping_ids(gqd.gene_wise_quantification._overlapping_alignments(
        sam, Gff3EntryMoc("chrom", 19, 23))) == [
            "myread:01", "myread:02", "myread:03", "myread:04", "myread:05"]
    # No overlap - very close downstream of the reads
    assert mapping_ids(gqd.gene_wise_quantification._overlapping_alignments(
        sam, Gff3EntryMoc("chrom", 20, 23))) == []


def test_overlapping_alignments_2():
    """Extraction of overlapping reads - with a non-default
        minimal overlap.
    """
    generate_bam_file(gqd.sam_content, gqd.sam_bam_prefix)
    gqd.gene_wise_quantification._min_overlap = 5
    sam = pysam.Samfile(gqd.sam_bam_prefix + ".bam")
    # 1 overlapping base in the 5' end of the reads => not enough
    assert mapping_ids(gqd.gene_wise_quantification._overlapping_alignments(
        sam, Gff3EntryMoc("chrom", 1, 10))) == []
    # 4 overlapping base in the 5' end of the reads => not enough
    assert mapping_ids(gqd.gene_wise_quantification._overlapping_alignments(
        sam, Gff3EntryMoc("chrom", 1, 13))) == []
    # 5 overlapping base in the 5' end of the reads => okay
    assert mapping_ids(gqd.gene_wise_quantification._overlapping_alignments(
        sam, Gff3EntryMoc("chrom", 1, 14))) == [
            "myread:01", "myread:02", "myread:03", "myread:04", "myread:05"]
    # 1 overlapping base in the 3' end of the reads => not enough
    assert mapping_ids(gqd.gene_wise_quantification._overlapping_alignments(
        sam, Gff3EntryMoc("chrom", 19, 23))) == []
    # 4 overlapping base in the 3' end of the reads => not enough
    assert mapping_ids(gqd.gene_wise_quantification._overlapping_alignments(
        sam, Gff3EntryMoc("chrom", 16, 23))) == []
    # 5 overlapping base in the 3' end of the reads => not enough
    assert mapping_ids(gqd.gene_wise_quantification._overlapping_alignments(
        sam, Gff3EntryMoc("chrom", 15, 23))) == [
            "myread:01", "myread:02", "myread:03", "myread:04", "myread:05"]
    

def mapping_ids(mappings):
    return [mapping.qname for mapping in mappings]


def generate_bam_file(sam_content, file_prefix):
    sam_file = "{}.sam".format(file_prefix)
    bam_file = "{}.bam".format(file_prefix)
    sam_fh = open(sam_file, "w")
    sam_fh.write(sam_content)
    sam_fh.close()
    pysam.view(
        "-Sb", "-o{}".format(bam_file), sam_file, catch_stdout=False)
    pysam.index(bam_file)
    #bam = pysam.Samfile(bam_file)
    #return bam


class Gff3EntryMoc(object):

    def __init__(self, seq_id, start, end):
        self.seq_id = seq_id
        self.start = start
        self.end = end
