import sys
from functools import reduce
sys.path.append(".")
from libs.fasta import FastaParser
from libs.sam import SamParser

class ReadMapperStats(object):

    def __init__(self):
        self.fasta_parser = FastaParser()
        self.sam_parser = SamParser()

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

    def count_mappings(self, read_file_names, read_mapping_result_paths):
        for read_file_name, read_mapping_result_path in zip(
            read_file_names, read_mapping_result_paths):
            self._count_mappings(open(read_mapping_result_path))
            
    def _count_mappings(self, read_mapping_result_fh):
        ref_seqs_and_mappings = {}
        ref_seqs_and_mapped_reads = {}
        for ref_seq in self.sam_parser.reference_sequences(
            read_mapping_result_fh):
            ref_seqs_and_mappings[ref_seq] = 0
            ref_seqs_and_mapped_reads[ref_seq] = 0
        for entry in self.sam_parser.entries(read_mapping_result_fh):
            try:
                ref_seqs_and_mappings[entry.reference] += 1
                ref_seqs_and_mapped_reads[
                    entry.reference] += 1.0/float(entry.number_of_hits_as_int)
            except:
                sys.stderr.write(
                    "SAM entry with unspecified reference found! Stoping\n")
                sys.exit(2)
        return(ref_seqs_and_mappings, ref_seqs_and_mapped_reads)

    def count_unmapped_reads(self):
        pass

    def _count_fasta_fh_entries(self, fasta_fh):
        # A memory saving approach to sum the number of entries
        return(reduce(lambda x, y: x+1, 
                      self.fasta_parser.entries(fasta_fh), 0))

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
        
