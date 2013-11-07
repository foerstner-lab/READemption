import sys
from trapllib.fasta import FastaParser
from trapllib.polyaclipper import PolyAClipper

class ReadClipper(object):

    def __init__(self):
        self.fasta_parser = FastaParser()
        self.poly_a_clipper = PolyAClipper()

    def clip(self, input_paths, output_pathes):
        for input_path, output_path in zip(
            input_paths, output_pathes):
            with open(input_path) as input_fh:
                  with open(output_path, "w") as output_fh:
                      self._clip_entries_in_file(input_fh, output_fh)

    def _compare_input_paths(self):
        if len(input_paths) != len(output_pathes):
            raise Exception("Number of input != number of output files.")
    
    def _clip_entries_in_file(self, input_fh, output_fh):
        for header, seq in self.fasta_parser.entries(input_fh):
            clipped_seq = self.poly_a_clipper.clip_before_long_poly_a_strech(
                seq)
            clipped_seq = self.poly_a_clipper.remove_3_prime_a(clipped_seq)
            output_fh.write(">%s\n%s\n" % (header, clipped_seq))
