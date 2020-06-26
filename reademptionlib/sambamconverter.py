import os
import pysam


class SamToBamConverter(object):
    def __init__(self):
        self._unsorted_appendix = ".tmp_unsorted"

    def sam_to_bam(self, sam_path, bam_path_prefix):
        if self._sam_file_is_empty(sam_path) is True:
            # pysam will generate an error if an emtpy SAM file will
            # be converted. Due to this an empty bam file with the
            # same header information will be generated from scratch
            self._generate_empty_bam_file(sam_path, bam_path_prefix)
            # Remove SAM file
            os.remove(sam_path)
            return
        temp_unsorted_bam_path = self._temp_unsorted_bam_path(bam_path_prefix)
        # Generate unsorted BAM file
        pysam.samtools.view(
            "-b",
            "-o{}".format(temp_unsorted_bam_path),
            sam_path,
            catch_stdout=False,
        )
        # Generate sorted BAM file
        pysam.sort(temp_unsorted_bam_path, "-o", bam_path_prefix + ".bam")
        # Generate index for BAM file
        pysam.index("%s.bam" % bam_path_prefix)
        # Remove unsorted BAM file
        os.remove(temp_unsorted_bam_path)
        # Remove SAM file
        os.remove(sam_path)

    def bam_to_sam(self, bam_path, sam_path):
        pysam.view("-ho{}".format(sam_path), bam_path, catch_stdout=False)

    def _temp_unsorted_bam_path(self, bam_path_prefix):
        return "%s%s.bam" % (bam_path_prefix, self._unsorted_appendix)

    def _sam_file_is_empty(self, sam_path):
        # Check if there is any line that is not a header line
        # (i.e. which is not starting with @)
        for line in open(sam_path):
            if line.startswith("@") is False:
                return False
        return True

    def _generate_empty_bam_file(self, sam_path, bam_path_prefix):
        samfile = pysam.Samfile(sam_path, "r")
        bamfile = pysam.Samfile(
            "%s.bam" % bam_path_prefix, "wb", header=samfile.header
        )
        bamfile.close()
        samfile.close()
        pysam.index("%s.bam" % bam_path_prefix)
