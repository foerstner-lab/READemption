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
        call("%s view -Sb -o %s %s" % 
             (self._samtools_bin, temp_unsorted_bam_file_path, sam_file_path), 
             shell=True)
        # Generate sorted BAM file
        call("%s sort %s %s" % 
             (self._samtools_bin, temp_unsorted_bam_file_path, 
              bam_file_path_prefix), shell=True)
        # Remove unsoreted BAM file
        os.remove(temp_unsorted_bam_file_path)
        # Remove SAM file
        #os.remove(sam_file_path)

    def _temp_unsorted_bam_file_path(self, bam_file_path_prefix):
        return("%s%s.bam" % (bam_file_path_prefix, self._unsorted_appendix))


class BamToSamStreamer(object):
    
    def __init__(self, samtools_bin):
        self._samtools_bin = samtools_bin

    def bam_to_sam_stream(self, bam_file):
        return(Popen("%s view -h %s" % (self._samtools_bin, bam_file), 
                  stdout=PIPE, shell=True).stdout)
