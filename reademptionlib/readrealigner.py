import sys
from reademptionlib.lack import Lack

class ReadRealigner(object):
    """An abstraction layer for different short read realigners."""

    def __init__(self, lack_bin, show_progress):
        self.lack = Lack(lack_bin, show_progress=show_progress)
    
    def run_alignment(self, query_fasta_path, query_sam_path, ref_seq_paths,
                      output_sam_path, nomatch_path, threads, accuracy):
        self.lack.realign_reads(query_fasta_path, 
                                query_sam_path, 
                                ref_seq_paths, 
                                output_sam_path,
                                nonmatch_path=nomatch_path, threads=threads,
                                accuracy=accuracy)
