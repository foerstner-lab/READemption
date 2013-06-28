import gzip
import bz2
from collections import defaultdict
from libs.fasta import FastaParser
from libs.polyaclipper import PolyAClipper

class ReadProcessor(object):

    def __init__(self, poly_a_clipping=False,  min_read_length=12):
        self._poly_a_clipping = poly_a_clipping
        self._min_read_length = min_read_length
        self.fasta_parser = FastaParser()
        self.poly_a_clipper = PolyAClipper()

    def process(self, input_path, output_path):
        self._stats = defaultdict(int)
        self._stats["total_no_of_reads"]
        self._stats["polya_removed"]
        self._stats["single_a_removed"]
        self._stats["unmodified"]
        self._stats["too_short"]
        self._stats["long_enough"]
        self._stats["read_length_before_processing_and_freq"] = defaultdict(int)
        self._stats["read_length_after_processing_and_freq"] = defaultdict(int)
        with gzip.open(output_path, "wb") as output_fh:
            input_fh = self._input_fh(input_path)
            self._process(input_fh, output_fh)
        return self._stats

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

    def _process(self, input_fh, output_fh):
        for header, seq in self.fasta_parser.entries(input_fh):
            self._stats["total_no_of_reads"] += 1
            if self._poly_a_clipping:
                clipped_seq = self.poly_a_clipper.clip_poly_a_strech(seq)
                clipped_seq = self.poly_a_clipper.remove_3_prime_a(clipped_seq)
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
