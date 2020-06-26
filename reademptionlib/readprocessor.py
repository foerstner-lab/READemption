import gzip
import bz2
from collections import defaultdict
from reademptionlib.polyaclipper import PolyAClipper
from Bio import SeqIO
from Bio.Seq import Seq


class ReadProcessor(object):
    def __init__(
        self,
        poly_a_clipping=False,
        min_read_length=12,
        paired_end=False,
        fastq=False,
        min_phred_score=None,
        adapter=None,
        reverse_complement=False,
    ):
        self._poly_a_clipping = poly_a_clipping
        self._min_read_length = min_read_length
        self._paired_end = paired_end
        self._fastq = fastq
        self._min_phred_score = min_phred_score
        self._adapter = adapter
        self._poly_a_clipper = PolyAClipper()
        self._reverse_complement = reverse_complement

    def process_single_end(self, input_path, output_path):
        self._init_stat_dict()
        with gzip.open(output_path, "wb") as output_fh:
            input_fh = self._input_fh(input_path)
            self._process_single_end(input_fh, output_fh)
        return self._stats

    def process_paired_end(self, input_path_pair, output_path_pair):
        self._init_stat_dict()
        with gzip.open(output_path_pair[0], "wb") as output_p1_fh, gzip.open(
            output_path_pair[1], "wb"
        ) as output_p2_fh:
            input_p1_fh = self._input_fh(input_path_pair[0])
            input_p2_fh = self._input_fh(input_path_pair[1])
            self._process_paired_end(
                input_p1_fh, input_p2_fh, output_p1_fh, output_p2_fh
            )
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
            return gzip.open(input_path, "rt")
        elif input_path.endswith(".bz2"):
            return bz2.open(input_path, "rt")
        return open(input_path)

    def _trim_by_quality(self, seq, qualities):
        good_nucl = []
        for nucl, qual in zip(seq, qualities):
            if qual < self._min_phred_score:
                break
            good_nucl.append(nucl)
        return "".join(good_nucl)

    def _clip_adapter(self, seq):
        adapter_start_pos = seq.find(self._adapter)
        if adapter_start_pos == -1:
            return seq
        else:
            return seq[:adapter_start_pos]

    def _process_single_end(self, input_fh, output_fh):
        for header, seq, qualities in self._parse_sequences(input_fh):
            raw_seq_len = len(seq)
            self._stats["total_no_of_reads"] += 1
            if self._fastq and not self._min_phred_score is None:
                seq = self._trim_by_quality(seq, qualities)
            if self._reverse_complement:
                seq = Seq(seq)
                seq = str(seq.reverse_complement())
            if not self._adapter is None:
                seq = self._clip_adapter(seq)
            if self._poly_a_clipping:
                seq = self._poly_a_clipper.clip_poly_a_strech(seq)
                seq = self._poly_a_clipper.remove_3_prime_a(seq)
            clipped_seq_len = len(seq)
            if clipped_seq_len == raw_seq_len - 1:
                self._stats["single_a_removed"] += 1
            elif clipped_seq_len < raw_seq_len - 1:
                self._stats["polya_removed"] += 1
            else:
                self._stats["unmodified"] += 1
            if clipped_seq_len < self._min_read_length:
                self._stats["too_short"] += 1
                continue
            self._stats["long_enough"] += 1
            self._stats["read_length_before_processing_and_freq"][
                raw_seq_len
            ] += 1
            self._stats["read_length_after_processing_and_freq"][
                clipped_seq_len
            ] += 1
            # Encoding to bytes is necessary due to saving via gzip
            output_fh.write(str.encode(">%s\n%s\n" % (header, seq)))

    def _parse_sequences(self, input_fh):
        if self._fastq:
            for seq_record in SeqIO.parse(input_fh, "fastq"):
                yield (
                    seq_record.description,
                    str(seq_record.seq),
                    seq_record.letter_annotations["phred_quality"],
                )
        else:
            for seq_record in SeqIO.parse(input_fh, "fasta"):
                yield (seq_record.description, str(seq_record.seq), None)

    def _process_paired_end(
        self, input_p1_fh, input_p2_fh, output_p1_fh, output_p2_fh
    ):
        for fasta_entry_p1, fasta_entry_p2 in zip(
            self._parse_sequences(input_p1_fh),
            self._parse_sequences(input_p2_fh,),
        ):
            header_p1 = fasta_entry_p1[0]
            header_p2 = fasta_entry_p2[0]
            seq_p1 = fasta_entry_p1[1]
            seq_p2 = fasta_entry_p2[1]
            qualities_p1 = fasta_entry_p1[2]
            qualities_p2 = fasta_entry_p1[2]
            raw_seq_p1_len = len(seq_p1)
            raw_seq_p2_len = len(seq_p2)
            self._stats["total_no_of_reads"] += 1
            self._stats["unmodified"] += 1
            if self._fastq and self._min_phred_score is not None:
                seq_p1 = self._trim_by_quality(seq_p1, qualities_p1)
                seq_p2 = self._trim_by_quality(seq_p2, qualities_p2)
            if self._reverse_complement:
                seq_p1 = Seq(seq_p1)
                seq_p1 = str(seq_p1.reverse_complement())
                seq_p2 = Seq(seq_p2)
                seq_p2 = str(seq_p2.reverse_complement())
            if self._adapter is not None:
                seq_p1 = self._clip_adapter(seq_p1)
                seq_p2 = self._clip_adapter(seq_p2)
            if (
                raw_seq_p1_len < self._min_read_length
                or raw_seq_p2_len < self._min_read_length
            ):
                self._stats["too_short"] += 1
                continue
            self._stats["long_enough"] += 1
            self._stats["read_length_before_processing_and_freq"][
                raw_seq_p1_len
            ] += 1
            self._stats["read_length_after_processing_and_freq"][
                raw_seq_p1_len
            ] += 1
            self._stats["read_length_before_processing_and_freq"][
                raw_seq_p2_len
            ] += 1
            self._stats["read_length_after_processing_and_freq"][
                raw_seq_p2_len
            ] += 1
            # Encoding to bytes is necessary due to saving via gzip
            output_p1_fh.write(str.encode(">%s\n%s\n" % (header_p1, seq_p1)))
            output_p2_fh.write(str.encode(">%s\n%s\n" % (header_p2, seq_p2)))
