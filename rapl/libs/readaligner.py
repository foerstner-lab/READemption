import sys
sys.path.append(".")
from libs.segemehl import Segemehl

class ReadAligner(object):

    def __init__(self, segemehl_bin):
        self.segemehl = Segemehl(segemehl_bin)
    
    def build_index(self, ref_seq_paths, index_path):
        self.segemehl.build_index(ref_seq_paths, index_path)

    def run_alignment(self, read_path, index_path, ref_seq_path, output_path, 
                      nomatch_path, threads, accuracy, evalue, split):
        self.segemehl.align_reads(read_path, index_path, ref_seq_path,
                                  output_path, nonmatch_file=nomatch_path,
                                  threads=threads, accuracy=accuracy, 
                                  evalue=evalue, split=split)
