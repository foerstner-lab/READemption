import sys
sys.path.append(".")
from libs.fasta import FastaParser

class SeqSizeFilter(object):

    def __init__(self):
        self.fasta_parser = FastaParser()

    def filter(self, input_file_paths, output_long_enough_paths, 
               output_too_short_paths, min_seq_length):
        self._compare_input_paths(input_file_paths, output_long_enough_paths,
               output_too_short_paths)
        for input_file_path, output_long_enough_path, output_too_short_path in zip(
            input_file_paths, output_long_enough_paths, output_too_short_paths):
            input_fh = open(input_file_path)
            output_long_enough_fh = open(output_long_enough_path, "w")
            output_too_short_fh = open(output_too_short_path, "w")
            self._filter_entries_in_file(
                input_fh, output_long_enough_fh, output_too_short_fh, 
                min_seq_length)
            input_fh.close()
            output_long_enough_fh.close()
            output_too_short_fh.close()

    def _compare_input_paths(self, input_file_paths, output_long_enough_paths,
               output_too_short_paths):
        if (len(input_file_paths) != len(output_long_enough_paths) or
            len(input_file_paths) != len(output_too_short_paths)):
            raise Exception("Number of input != number of output files.")
    
    def _filter_entries_in_file(
        self, input_fh, output_long_enough_fh, output_too_short_fh, 
        min_seq_length):
        for header, seq in self.fasta_parser.entries(input_fh):
            output_fh = None
            if len(seq) >= min_seq_length:
                output_fh = output_long_enough_fh
            else:
                output_fh = output_too_short_fh
            output_fh.write(">%s\n%s\n" % (header, seq))
