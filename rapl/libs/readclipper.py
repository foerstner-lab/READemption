import sys
sys.path.append(".")
from libs.fasta import FastaParser
from libs.polyaclipper import PolyAClipper

class ReadClipper(object):

    def __init__(self):
        self.fasta_parser = FastaParser()
        self.poly_a_clipper = PolyAClipper()

    def clip(self, input_file_paths, output_file_pathes):
        for input_file_path, output_path in zip(
            input_file_paths, output_file_pathes):
            input_fh = open(input_file_path)
            output_fh = open(output_path, "w")
            self._clip_entries_in_file(input_fh, output_fh)
        input_fh.close()
        output_fh.close()

    def _compare_input_paths(self):
        if len(input_file_paths) != len(output_file_pathes):
            raise Exception("Number of input != number of output files.")
    
    def _clip_entries_in_file(self, input_fh, output_fh):
        for header, seq in self.fasta_parser.entries(input_fh):
            clipped_seq = self.poly_a_clipper.clip_before_long_poly_a_strech(
                seq)
            clipped_seq = self.poly_a_clipper.remove_3_prime_a(clipped_seq)
            output_fh.write(">%s\n%s\n" % (header, clipped_seq))
