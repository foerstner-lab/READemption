import gzip
import bz2
from collections import defaultdict
from reademptionlib.fasta import FastaParser
from reademptionlib.fastq import FastqParser
from reademptionlib.polyaclipper import PolyAClipper

class ReadProcessor(object):

    def __init__(self, poly_a_clipping=False,  min_read_length=12, 
                 paired_end=False, fastq=False):
        self._poly_a_clipping = poly_a_clipping
        self._min_read_length = min_read_length
        self._paired_end = paired_end
        self._fastq = fastq
        self._poly_a_clipper = PolyAClipper()

    def process_single_end(self, input_path, output_path):
        self._init_stat_dict()
        with gzip.open(output_path, "wb") as output_fh:
            input_fh = self._input_fh(input_path)
            self._process_single_end(input_fh, output_fh)
        return self._stats

    def process_paired_end(self, input_path_pair, output_path_pair):
        self._init_stat_dict()
        with gzip.open(output_path_pair[0], "wb") as output_p1_fh, \
                gzip.open(output_path_pair[1], "wb") as output_p2_fh:
            input_p1_fh = self._input_fh(input_path_pair[0])
            input_p2_fh = self._input_fh(input_path_pair[1])
            self._process_paired_end(
                input_p1_fh, input_p2_fh, output_p1_fh, output_p2_fh)
        return self._stats

    def _init_stat_dict(self):
        self._stats = defaultdict(int)
        self._stats["total_no_of_reads"]
        self._stats["polya_removed"]
        self._stats["single_a_removed"]
        self._stats["unmodified"]
        self._stats["too_short"]
        self._stats["long_enough"]
        self._stats["read_length_before_processing_and_freq"] = defaultdict(int)
        self._stats["read_length_after_processing_and_freq"] = defaultdict(int)        
    
    def _input_fh(self, input_path):
        """Return a file hande 

        Can deal with plain fasta files, gzipped fasta or bzipped2 fasta.
        """
        if input_path.endswith(".gz"):
            return self._line_interator(input_path, gzip.open)
        elif input_path.endswith(".bz2"):
            return self._line_interator(input_path, bz2.open)
        return open(input_path)

    def _line_interator(self, input_path, open_func):
        """Return a line iterator that decodes the bytes"""
        with open_func(input_path, "rb") as fh:
            for line in fh:
                yield(line.decode())

    def _process_single_end(self, input_fh, output_fh):
        for header, seq in self._seq_parser().entries(input_fh):
            self._stats["total_no_of_reads"] += 1
            if self._poly_a_clipping:
                clipped_seq = self._poly_a_clipper.clip_poly_a_strech(seq)
                clipped_seq = self._poly_a_clipper.remove_3_prime_a(clipped_seq)
            else:
                clipped_seq = seq
            if len(clipped_seq) == len(seq) - 1:
                self._stats["single_a_removed"] += 1
            elif len(clipped_seq) < len(seq) - 1:
                self._stats["polya_removed"] += 1
            else:
                self._stats["unmodified"] += 1
            clipped_seq_len = len(clipped_seq)
            if clipped_seq_len < self._min_read_length:
                self._stats["too_short"] += 1
                continue
            self._stats["long_enough"] += 1
            raw_seq_len = len(seq)
            self._stats["read_length_before_processing_and_freq"][
                raw_seq_len] += 1
            self._stats["read_length_after_processing_and_freq"][
                clipped_seq_len] += 1
            # Encoding to bytes is necessary due to saving via gzip
            output_fh.write(str.encode(">%s\n%s\n" % (header, clipped_seq)))

    def _process_paired_end(
        self, input_p1_fh, input_p2_fh, output_p1_fh, output_p2_fh):
        for fasta_entry_p1, fasta_entry_p2 in zip(
            self._seq_parser().entries(input_p1_fh), 
            self._seq_parser().entries(input_p2_fh)):
            header_p1 = fasta_entry_p1[0]
            header_p2 = fasta_entry_p2[0]
            seq_p1 = fasta_entry_p1[1]
            seq_p2 = fasta_entry_p2[1]
            self._stats["total_no_of_reads"] += 1
            self._stats["unmodified"] += 1
            seq_p1_len = len(seq_p1)
            seq_p2_len = len(seq_p2)
            if (seq_p1_len < self._min_read_length or 
                seq_p2_len < self._min_read_length):
                self._stats["too_short"] += 1
                continue
            self._stats["long_enough"] += 1
            self._stats["read_length_before_processing_and_freq"][
                seq_p1_len] += 1
            self._stats["read_length_after_processing_and_freq"][
                seq_p1_len] += 1
            self._stats["read_length_before_processing_and_freq"][
                seq_p2_len] += 1
            self._stats["read_length_after_processing_and_freq"][
                seq_p2_len] += 1
            # Encoding to bytes is necessary due to saving via gzip
            output_p1_fh.write(str.encode(">%s\n%s\n" % (header_p1, seq_p1)))
            output_p2_fh.write(str.encode(">%s\n%s\n" % (header_p2, seq_p2)))

    def _seq_parser(self):
        if self._fastq: return FastqParser()
        else: return FastaParser()

