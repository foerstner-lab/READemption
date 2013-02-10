from libs.fasta import FastaParser
from libs.polyaclipper import PolyAClipper

class ReadProcessor(object):

    def __init__(self, poly_clipping=True,  min_read_length=12):
        self.poly_clipping = poly_clipping
        self._min_read_length = min_read_length
        self.fasta_parser = FastaParser()
        self.poly_a_clipper = PolyAClipper()

    def process(self, input_path, output_path):
        self._stats = {
            "total_no_of_reads" : 0,
            "polya_removed" : 0,
            "single_a_removed" : 0,
            "unmodified" : 0,
            "too_short" : 0,
            "long_enough" : 0,
            "read_length_before_processing_and_freq" : {},
            "read_length_after_processing_and_freq" : {}}
        output_fh = open(output_path, "w")
        self._process(open(input_path), output_fh)
        output_fh.close()
        return(self._stats)

    def _process(self, input_fh, output_fh):
        for header, seq in self.fasta_parser.entries(input_fh):
            self._stats["total_no_of_reads"] += 1
            if self.poly_clipping:
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
            self._stats["read_length_before_processing_and_freq"].setdefault(
                raw_seq_len, 0)
            self._stats["read_length_before_processing_and_freq"][
                raw_seq_len] += 1
            self._stats["read_length_after_processing_and_freq"].setdefault(
                clipped_seq_len, 0)
            self._stats["read_length_after_processing_and_freq"][
                clipped_seq_len] += 1
            output_fh.write(">%s\n%s\n" % (header, clipped_seq))
