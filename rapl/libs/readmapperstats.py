import csv
import sys
from functools import reduce
sys.path.append(".")
from libs.fasta import FastaParser
from libs.sam import SamParser

class ReadMapperStats(object):

    def __init__(self, samtools_bin="samtools"):
        self.fasta_parser = FastaParser()
        self.sam_parser = SamParser(samtools_bin)

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

    def count_mappings(self, read_file_names, read_mapping_result_bam_paths):
        # self.no_of_mappings, self.no_of_mapped_reads and
        # self.no_of_uniquely_mapped_reads are dictionaries of
        # dictionaries: Read file name -> Reference seq -> counting
        self.no_of_mappings = {} 
        self.no_of_mapped_reads = {}
        self.no_of_uniquely_mapped_reads = {}
        for read_file_name, read_mapping_result_bam_path in zip(
            read_file_names, read_mapping_result_bam_paths):
            ref_seq_ids = self.sam_parser.ref_seq_ids_and_lengths_bam(
                  read_mapping_result_bam_path).keys()
            (self.no_of_mappings[read_file_name], 
             self.no_of_mapped_reads[read_file_name], 
             self.no_of_uniquely_mapped_reads[read_file_name]) = (
                self.sam_parser.mapping_countings_bam(
                     read_mapping_result_bam_path, ref_seq_ids))

    def count_unmapped_reads(self, read_file_names, unmapped_read_paths):
        self.no_of_unmapped_reads = {}
        self._count_fasta_entry_set(
            read_file_names, unmapped_read_paths, self.no_of_unmapped_reads)

    def _count_fasta_fh_entries(self, fasta_fh):
        # A memory saving approach to sum the number of entries
        return(reduce(lambda x, y: x+1, 
                      self.fasta_parser.entries(fasta_fh), 0))

    def write_stats_to_file(
        self, read_file_names, ref_ids_to_file_name, output_file_path):
        self._write_stats_to_fh(
            read_file_names, ref_ids_to_file_name, open(output_file_path, "w"))

    def _write_stats_to_fh(
        self, read_file_names, ref_ids_to_file_name, output_fh):
        output_fh.write(self._head_line(read_file_names) + "\n")
        for description, value_dict in [
            ("Number of raw reads", self.raw_read_countings),
            ("Reads long enough after clipping", 
             self.long_enough_clipped_reads),
            ("Reads too short after clipping", 
             self.too_small_clipped_reads)]: 
            output_fh.write(self._value_line(
                    description, value_dict, read_file_names) + "\n")
        for description, value_dict_of_dicts in [
            ("Total number of mapped reads", self.no_of_mapped_reads),
            ("Total number of uniquely mapped reads", 
             self.no_of_uniquely_mapped_reads),
            ("Total number of mappings", self.no_of_mappings)]:
            output_fh.write(self._dict_value_sum_line(
                    description, value_dict_of_dicts, read_file_names) + "\n")
        output_fh.write(self._value_line(
                "Number of unmappped reads", 
                self.no_of_unmapped_reads, read_file_names) + "\n")
        ref_seq_headers = sorted(
            list(self.no_of_mapped_reads.items())[0][1].keys())
        for ref_seq_header in ref_seq_headers:
            output_fh.write(
                self._dict_value_per_ref_genome_line(
                    "Number of mapped reads in %s" % 
                    ref_ids_to_file_name[ref_seq_header],
                    self.no_of_mapped_reads, read_file_names, ref_seq_header)
                + "\n")
        for ref_seq_header in ref_seq_headers:
            output_fh.write(
                self._dict_value_per_ref_genome_line(
                    "Number of uniquely mapped reads in %s" % 
                    ref_ids_to_file_name[ref_seq_header],
                    self.no_of_uniquely_mapped_reads, read_file_names, 
                    ref_seq_header) + "\n")
        for ref_seq_header in ref_seq_headers:
            output_fh.write(
                self._dict_value_per_ref_genome_line(
                    "Number of mapping in %s" % 
                    ref_ids_to_file_name[ref_seq_header],
                    self.no_of_mappings, read_file_names, ref_seq_header)
                + "\n")

    def _head_line(self, read_file_names):
        return("\t" + "\t".join(read_file_names))

    def _value_line(self, description, value_dict, read_file_names):
        return(description + "\t" + "\t".join(
                [str(value_dict[read_file_name])
                 for read_file_name in read_file_names]))
        
    def _dict_value_sum_line(
        self, description, value_dict_of_dicts, read_file_names):
        return(description + "\t" + "\t".join(
                [str(sum(value_dict_of_dicts[read_file_name].values()))
                 for read_file_name in read_file_names]))

    def _dict_value_per_ref_genome_line(
        self, description, value_dict_of_dicts, read_file_names, ref_seq_header):
        return(description + "\t" + "\t".join(
                [str(value_dict_of_dicts[read_file_name][ref_seq_header])
                 for read_file_name in read_file_names]))


class ReadMapperStatsReader(object):
    
    def read_stat_file(self, stat_file_path):
        return(self._read_stat_file(open(stat_file_path)))

    def _read_stat_file(self, stat_fh):
        stats = {}
        libs = stat_fh.readline()[:-1].split("\t")[1:]
        total_num_of_mapped_read = None
        for row in csv.reader(stat_fh, delimiter="\t"):
            if row[0].startswith("Total number of mapped reads"):
                total_num_of_mapped_read = [float(count) for count in row[1:]]
                break

        for lib, counting in zip(libs, total_num_of_mapped_read):
            stats[lib] = {"total_number_of_mapped_reads" : counting}
        return(stats)

    def min_read_countings(self, stat_file_path):
        read_mapping_stats = self.read_stat_file(stat_file_path)
        return(min([lib_features["total_number_of_mapped_reads"] 
                    for lib_features in read_mapping_stats.values()]))
