import sys
from functools import reduce
sys.path.append(".")
from libs.fasta import FastaParser

class ReadMapperStats(object):

    def count_raw_reads(self, read_file_names, read_file_paths):
        self.raw_read_countings = {}
        self._count_fasta_entry_set(
            read_file_names, read_file_paths, self.raw_read_countings)

    def count_too_small_clipped_reads(self, read_file_names, read_file_paths):
        self.too_small_clipped_reads = {}
        self._count_fasta_entry_set(
            read_file_names, read_file_paths, self.too_small_clipped_reads)

    def count_long_enough_clipped_reads(
        self, read_file_names, read_file_paths):
        self.long_enough_clipped_reads = {}
        self._count_fasta_entry_set(
            read_file_names, read_file_paths, self.long_enough_clipped_reads)

    def _count_fasta_entry_set(self, file_names, file_paths, saving_dict):
        for file_name, file_path in zip(file_names, file_paths):
            saving_dict[file_name] = self._count_fasta_entries(file_path)

    def _count_fasta_entries(self, fasta_path):
        return(self._count_fasta_fh_entries(open(fasta_path)))

    def _count_fasta_fh_entries(self, fasta_fh):
        fasta_parser = FastaParser()
        # A memory saving approach to sum the number of entries
        return(reduce(lambda x, y: x+1, fasta_parser.entries(fasta_fh), 0))

    def write_stats_to_file(self, read_file_names, output_file_path):
        self._write_stats_to_fh(read_file_names, open(output_file_path, "w"))

    def _write_stats_to_fh(self, read_file_names, output_fh):
        output_fh.write(self._head_line(read_file_names) + "\n")
        for description, value_dict in [
            ("Number of raw reads", self.raw_read_countings),
            ("Reads long enough after clipping", 
             self.long_enough_clipped_reads),
            ("Reads too short after clipping", 
             self.too_small_clipped_reads)]:
            output_fh.write(self._value_line(
                    description, value_dict, read_file_names) + "\n")

    def _head_line(self, read_file_names):
        return("\t" + "\t".join(read_file_names))

    def _value_line(self, description, value_dict, read_file_names):
        return(description + "\t" + "\t".join(
                [str(value_dict[read_file_name])
                 for read_file_name in read_file_names]))
        
