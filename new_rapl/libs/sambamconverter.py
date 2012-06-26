from subprocess import call, Popen, PIPE
import os

class SamToBamConverter(object):
    
    def __init__(self, samtools_bin):
        self._samtools_bin = samtools_bin
        self._unsorted_appendix = ".tmp_unsorted"
        
    def sam_to_bam(self, sam_file_path, bam_file_path_prefix):
        temp_unsorted_bam_file_path = self._temp_unsorted_bam_file_path(
            bam_file_path_prefix)
        # Generate unsorted BAM file
        call([self._samtools_bin, "view", "-Sb", "-o", 
              temp_unsorted_bam_file_path, sam_file_path])
        # Generate sorted BAM file
        call([self._samtools_bin, "sort", temp_unsorted_bam_file_path, 
              bam_file_path_prefix])

        # Generate index for BAM file
        call("%s index %s.bam" % 
             (self._samtools_bin, bam_file_path_prefix), shell=True)
        # Remove unsorted BAM file
        os.remove(temp_unsorted_bam_file_path)
        # Remove SAM file
        os.remove(sam_file_path)

    def _temp_unsorted_bam_file_path(self, bam_file_path_prefix):
        return("%s%s.bam" % (bam_file_path_prefix, self._unsorted_appendix))
