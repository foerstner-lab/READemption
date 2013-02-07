import csv
import sys
from functools import reduce
sys.path.append(".")
from libs.fasta import FastaParser
import pysam

class ReadMapperStats(object):

    def __init__(self):
        self.fasta_parser = FastaParser()

    def count(self, read_mapping_result_bam_path, unmapped_reads_path):
        self._stats = {}
        self._count_mapped_reads_and_mappings(read_mapping_result_bam_path)
        self._count_unmapped_reads(unmapped_reads_path)
        return(self._stats)

    def _count_unmapped_reads(self, unmapped_read_paths):
        fasta_fh = open(unmapped_read_paths)
        self._stats["no_of_unmapped_reads"] = self._count_fasta_entries(
            fasta_fh)
        fasta_fh.close()

    def _count_fasta_entries(self, fasta_fh):
        return(reduce(lambda x, y: x + 1,
                      self.fasta_parser.entries(fasta_fh), 0))

    def _count_mapped_reads_and_mappings(self, read_mapping_result_bam_path):
        bam = pysam.Samfile(read_mapping_result_bam_path)
        stats_per_ref = {}
        no_of_hits_per_read_freq = {}
        for ref_id in bam.references:
            self._init_counting_dict(stats_per_ref, ref_id)
        for entry in bam.fetch():
            ref_id = bam.getrname(entry.tid)
            try:
                self._count_mapping(
                    entry, ref_id, stats_per_ref, no_of_hits_per_read_freq)
            except KeyError:
                sys.stderr.write(
                    "SAM entry with unspecified reference found! Stoping\n")
                sys.exit(2)

        self._stats["countings_per_reference"] = stats_per_ref
        self._stats["countings_total"] = self._sum_countings(stats_per_ref)
        self._stats["no_of_hits_per_read_and_freq"] = self._calc_down_to_read(
            no_of_hits_per_read_freq)

    def _sum_countings(self, stats_per_ref):
        total_stats = {}
        for ref_id, stats in stats_per_ref.items():
            for attribute, value in stats.items():
                total_stats.setdefault(attribute, 0)
                total_stats[attribute] += value
        return(total_stats)

    def _calc_down_to_read(self, no_of_hits_per_read_freq):
        """As the frequencies were determined via the mappings we need
        to normalized each frequency value down to the read by
        dividing the frequencig by the number of hits per read.
        """
        return(dict((no_of_hits_per_read, freq/no_of_hits_per_read)
                    for no_of_hits_per_read, freq in
                    no_of_hits_per_read_freq.items()))

    def _init_counting_dict(self, stats_per_ref, ref_id):
        stats_per_ref[ref_id] = {}
        stats_per_ref[ref_id]["no_of_mappings"] = 0
        stats_per_ref[ref_id]["no_of_mapped_reads"] = 0
        stats_per_ref[ref_id]["no_of_uniquely_mapped_reads"] = 0

    def _count_mapping(
            self, entry, ref_id, stats_per_ref, no_of_hits_per_read_freq):
        no_of_hits = dict(entry.tags)["NH"]
        no_of_hits_per_read_freq.setdefault(no_of_hits, 0)
        no_of_hits_per_read_freq[no_of_hits] += 1
        stats_per_ref[ref_id]["no_of_mappings"] += 1
        stats_per_ref[
            ref_id]["no_of_mapped_reads"] += 1.0/float(no_of_hits)
        if no_of_hits == 1:
            stats_per_ref[ref_id]["no_of_uniquely_mapped_reads"] += 1

class ReadMapperStatsReader(object):

    def read_mapping_stat_file(self, stat_file_path):
        return(self._read_stat_file(open(stat_file_path)))

    def _read_mapping_stat_file(self, stat_fh):
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
