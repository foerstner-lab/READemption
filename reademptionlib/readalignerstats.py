import csv
import sys
from collections import defaultdict
from functools import reduce
from reademptionlib.fasta import FastaParser
import pysam

class ReadAlignerStats(object):

    def __init__(self):
        self.fasta_parser = FastaParser()

    def count(self, read_alignment_result_bam_path, unaligned_reads_path):
        self._stats = {}
        self._count_aligned_reads_and_alignments(read_alignment_result_bam_path)
        self._count_unaligned_reads(unaligned_reads_path)
        return self._stats

    def _count_unaligned_reads(self, unaligned_read_paths):
        
        with open(unaligned_read_paths) as fasta_fh:
            self._stats["stats_total"][
                "no_of_unaligned_reads"] = self._count_fasta_entries(fasta_fh)

    def _count_fasta_entries(self, fasta_fh):
        return reduce(lambda x, y: x + 1,
                      self.fasta_parser.entries(fasta_fh), 0)

    def _count_aligned_reads_and_alignments(
            self, read_alignment_result_bam_path):
        bam = pysam.Samfile(read_alignment_result_bam_path)
        stats_per_ref = defaultdict(dict)
        no_of_hits_per_read_freq = {}
        for ref_id in bam.references:
            self._init_counting_dict(stats_per_ref, ref_id)
        for entry in bam.fetch():
            ref_id = bam.getrname(entry.tid)
            try:
                self._count_alignment(
                    entry, ref_id, stats_per_ref, no_of_hits_per_read_freq)
            except KeyError:
                sys.stderr.write(
                    "SAM entry with unspecified reference found! Stoping\n")
                sys.exit(2)
        self._stats["stats_per_reference"] = stats_per_ref
        for ref_id, stats in stats_per_ref.items():
            stats_per_ref[ref_id][
                "no_of_hits_per_read_and_freqs"] = self._calc_down_to_read(
                    stats_per_ref[ref_id]["no_of_hits_per_read_and_freqs"])
        self._stats["stats_total"] = self._sum_countings(stats_per_ref)

    def _sum_countings(self, stats_per_ref):
        total_stats = {}
        for ref_id, stats in stats_per_ref.items():
            for attribute, value in stats.items():
                if type(value) is int or type(value) is float:
                    total_stats.setdefault(attribute, 0)
                    total_stats[attribute] += value
                elif type(value) is dict:
                    total_stats.setdefault(attribute, {})
                    for value_int, freq in value.items():
                        total_stats[attribute].setdefault(value_int, 0)
                        total_stats[attribute][value_int] += freq
        return total_stats

    def _calc_down_to_read(self, no_of_hits_per_read_freq):
        """As the frequencies were determined via the alignments we need
        to normalized each frequency value down to the read by
        dividing the frequencig by the number of hits per read.
        """
        return dict((no_of_hits_per_read, freq/no_of_hits_per_read)
                    for no_of_hits_per_read, freq in
                    no_of_hits_per_read_freq.items())

    def _init_counting_dict(self, stats_per_ref, ref_id):
        stats_per_ref[ref_id] = defaultdict(float)
        stats_per_ref[ref_id]["no_of_alignments"]
        stats_per_ref[ref_id]["no_of_aligned_reads"]
        stats_per_ref[ref_id]["no_of_split_alignments"]
        stats_per_ref[ref_id]["no_of_uniquely_aligned_reads"]
        stats_per_ref[ref_id][
            "alignment_length_and_freqs"] = defaultdict(int)
        stats_per_ref[ref_id][
            "no_of_hits_per_read_and_freqs"] = defaultdict(int)

    def _count_alignment(self, entry, ref_id, stats_per_ref,
                         no_of_hits_per_read_freq):
        entry_tags_dict = dict(entry.tags)
        no_of_hits = entry_tags_dict["NH"]
        # Consider split reads
        no_of_splits = float(entry_tags_dict.get("XL", 1))
        stats_per_ref[ref_id]["no_of_hits_per_read_and_freqs"][
            no_of_hits] += 1
        if "XL" in entry_tags_dict:
            stats_per_ref[ref_id]["no_of_split_alignments"] += 1.0/no_of_splits
        stats_per_ref[ref_id]["no_of_alignments"] += 1.0/no_of_splits
        stats_per_ref[
            ref_id]["no_of_aligned_reads"] += 1.0/(
            float(no_of_hits) * no_of_splits)
        if no_of_hits == 1:
            stats_per_ref[ref_id][
                "no_of_uniquely_aligned_reads"] += 1.0/no_of_splits
        stats_per_ref[ref_id][
            "alignment_length_and_freqs"][entry.alen] += 1
