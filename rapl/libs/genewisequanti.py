import csv
import os.path
from libs.gff3 import Gff3Parser
import pysam

class GeneWiseQuantification(object):

    def __init__(self, min_overlap=1, norm_by_alignment_freq=True,
                 norm_by_overlap_freq=True):
        """
        - normalize_by_alignment: consider that some reads are aligned at
          more than one location and only count fractions
        - normalize_by_overlapping_genes: consider that some alignment
          overlap with more than on gene

        """
        self._min_overlap = min_overlap
        self._norm_by_alignment_freq = norm_by_alignment_freq
        self._norm_by_overlap_freq = norm_by_overlap_freq

    def calc_overlaps_per_alignment(self, read_alignment_path,
                                  annotation_paths):
        """Calculate for each alignment the number of genes it
        overlaps. This has to be done globally i.e. for all annotation
        files combined in one dictionary.
        """
        gff3_parser = Gff3Parser()
        self.alignments_and_no_of_overlaps = {}
        for annotation_path in annotation_paths:
            sam = pysam.Samfile(read_alignment_path)
            for entry in gff3_parser.entries(open(annotation_path)):
                for alignment in self._overlapping_alignments(sam, entry):
                    alignment_id = self._alignment_id(alignment)
                    self.alignments_and_no_of_overlaps.setdefault(alignment_id, 0)
                    self.alignments_and_no_of_overlaps[alignment_id] += 1

    def quantify(self, read_alignment_path, annotation_path, output_path):
        self._quantify(read_alignment_path, annotation_path, output_path,
                      self._fraction_calc_method())

    def _quantify(self, read_alignment_path, annotation_path, output_path,
                  fraction_calc_method):
        sam = pysam.Samfile(read_alignment_path)
        gff3_parser = Gff3Parser()
        output_fh = open(output_path, "w")
        output_fh.write("#" + "\t".join([""] * 9) + "sense\tantisense\n")
        for entry in gff3_parser.entries(open(annotation_path)):
            sum_sense = 0
            sum_antisense = 0
            for alignment in self._overlapping_alignments(sam, entry):
                alignment_id = self._alignment_id(alignment)
                fraction = fraction_calc_method(alignment)
                if self._same_strand(entry, alignment):
                    sum_sense += fraction
                else:
                    sum_antisense += fraction
            output_fh.write(str(entry) + "\t" + str(sum_sense) + "\t" +
                            str(sum_antisense) + "\n")

    def _same_strand(self, entry, alignment):
        assert entry.strand in ["+", "-"]
        if ((entry.strand == "+" and not alignment.is_reverse) or
            (entry.strand == "-" and alignment.is_reverse)):
            return(True)
        return(False)

    def _fraction_calc_method(self):
        if self._norm_by_alignment_freq and self._norm_by_overlap_freq:
           return(self._fraction_norm_by_alignment_and_overlap)
        elif self._norm_by_alignment_freq and not self._norm_by_overlap_freq:
            return(self._fraction_norm_by_alignment)
        elif not self._norm_by_alignment_freq and self._norm_by_overlap_freq:
            return(self._fraction_norm_by_overlap)
        return(self._fraction_calc_constant_one)

    def _fraction_calc_constant_one(self, alignment):
        return(1.0)

    def _fraction_norm_by_alignment_and_overlap(self, alignment):
        return(1.0 /
               float(self.alignments_and_no_of_overlaps[
                   self._alignment_id(alignment)]) /
                   float(dict(alignment.tags)["NH"])) # no. of alignments of read

    def _fraction_norm_by_alignment(self, alignment):
        return(1.0 / float(dict(alignment.tags)["NH"])) # no. of alignments of read

    def _fraction_norm_by_overlap(self, alignment):
        return(1.0 /
               float(self.alignments_and_no_of_overlaps[
                   self._alignment_id(alignment)]))

    def _overlapping_alignments(self, sam, entry):
        # The substraction of 1 from the start is necessary to perform
        # this correctly (checked in IGB, IGV and the unit testings).
        for alignment in sam.fetch(
                reference=entry.seq_id, start=entry.start-1, end=entry.end):
            if alignment.overlap(entry.start-1, entry.end) < self._min_overlap:
                continue
            yield(alignment)

    def _alignment_id(self, alignment):
        return(":".join([str(alignment.tid), alignment.qname, str(alignment.flag),
                         str(alignment.pos), str(alignment.aend)]))

    def _values_to_gene_key(self, seq_id, feature, start, end, strand):
        return("|".join(
                [str(val) for val in [seq_id, feature, start, end, strand]]))

class GeneWiseOverview(object):

    def create_overview(
            self, path_and_name_combos, read_file_names, overview_path):
        output_fh = open(overview_path, "w")
        # Write header
        output_fh.write("\t".join([""] * 10 + read_file_names) + "\n")
        self._add_to_overview(path_and_name_combos, "sense", 9, output_fh)
        self._add_to_overview(path_and_name_combos, "anti-sense", 10, output_fh)

    def _add_to_overview(self, path_and_name_combos, direction, column,
                         output_fh):
        gff3_parser = Gff3Parser()
        for annotation_path in sorted(path_and_name_combos.keys()):
            table_columns = []
            entries = []
            for entry in gff3_parser.entries(open(annotation_path)):
                entries.append(direction + "\t" + str(entry))
            table_columns.append(entries)
            for read_file, gene_quanti_path in path_and_name_combos[
                    annotation_path]:
                reader = csv.reader(open(gene_quanti_path), delimiter="\t")
                next(reader) # skip first line
                table_columns.append([row[column] for row in reader])
            # Generate a table by rotating the column list
            table = zip(*table_columns)
            for row in table:
                output_fh.write("\t".join(row) + "\n")
