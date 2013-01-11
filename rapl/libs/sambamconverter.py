import os
import pysam

class SamToBamConverter(object):

    def __init__(self):
        self._unsorted_appendix = ".tmp_unsorted"

    def sam_to_bam(self, sam_file_path, bam_file_path_prefix):
        temp_unsorted_bam_file_path = self._temp_unsorted_bam_file_path(
            bam_file_path_prefix)
        # Generate unsorted BAM file
        pysam.view("-Sb", "-o%s" % temp_unsorted_bam_file_path, sam_file_path)
        # Generate sorted BAM file
        pysam.sort(temp_unsorted_bam_file_path, bam_file_path_prefix)
        # Generate index for BAM file
        pysam.index("%s.bam" % bam_file_path_prefix)
        # Remove unsorted BAM file
        os.remove(temp_unsorted_bam_file_path)
        # Remove SAM file
        os.remove(sam_file_path)

    def _temp_unsorted_bam_file_path(self, bam_file_path_prefix):
        return("%s%s.bam" % (bam_file_path_prefix, self._unsorted_appendix))
