import csv
import os.path
from libs.gff3 import Gff3Parser
import pysam

class GeneWiseQuantification(object):

    def __init__(self, min_overlap=1, norm_by_mapping_freq=True,
                 norm_by_overlap_freq=True):
        """
        - normalize_by_mapping: consider that some reads are mapped at
          more than one location and only count fractions
        - normalize_by_overlapping_genes: consider that some mapping
          overlap with more than on gene

        """
        self._min_overlap = min_overlap
        self._norm_by_mapping_freq = norm_by_mapping_freq
        self._norm_by_overlap_freq = norm_by_overlap_freq

    def calc_overlaps_per_mapping(self, read_mapping_path,
                                  annotation_file_paths):
        """Calculate for each mapping the number of genes it
        overlaps. This has to be done for globally i.e. for all
        annotation files combined in one dictionary.
        """
        gff3_parser = Gff3Parser()
        self.mappings_and_no_of_overlaps = {}
        for annotation_file_path in annotation_file_paths:
            sam = pysam.Samfile(read_mapping_path)
            for entry in gff3_parser.entries(open(annotation_file_path)):
                for mapping in self._overlapping_mappings(sam, entry):
                    mapping_id = self._mapping_id(mapping)
                    self.mappings_and_no_of_overlaps.setdefault(mapping_id, 0)
                    self.mappings_and_no_of_overlaps[mapping_id] += 1

    def quantify(self, read_mapping_path, annotation_file_path, output_path):
        self._quantify(read_mapping_path, annotation_file_path, output_path,
                      self._fraction_calc_method())

    def _quantify(self, read_mapping_path, annotation_file_path, output_path,
                  fraction_calc_method):
        sam = pysam.Samfile(read_mapping_path)
        gff3_parser = Gff3Parser()
        sum_sense = 0
        sum_antisense = 0
        output_fh = open(output_path, "w")
        output_fh.write("#" + "\t".join([""] * 9) + "sense\tantisense\n")
        for entry in gff3_parser.entries(open(annotation_file_path)):
            for mapping in self._overlapping_mappings(sam, entry):
                mapping_id = self._mapping_id(mapping)
                fraction = fraction_calc_method(mapping)
                if self._same_strand(entry, mapping):
                    sum_sense += fraction
                else:
                    sum_antisense += fraction
            output_fh.write(str(entry) + "\t" + str(sum_sense) + "\t" +
                            str(sum_antisense) + "\n")

    def _same_strand(self, entry, mapping):
        assert entry.strand in ["+", "-"]
        if ((entry.strand == "+" and not mapping.is_reverse) or
            (entry.strand == "-" and mapping.is_reverse)):
            return(True)
        return(False)

    def _fraction_calc_method(self):
        if self._norm_by_mapping_freq and self._norm_by_overlap_freq:
           return(self._fraction_norm_by_mapping_and_overlap)
        elif self._norm_by_mapping_freq and not self._norm_by_overlap_freq:
            return(self._fraction_norm_by_mapping)
        elif not self._norm_by_mapping_freq and self._norm_by_overlap_freq:
            return(self._fraction_norm_by_overlap)
        return(self._fraction_calc_constant_one)

    def _fraction_calc_constant_one(self, mapping):
        return(1.0)

    def _fraction_norm_by_mapping_and_overlap(self, mapping):
        return(1.0 /
               float(self.mappings_and_no_of_overlaps[
                   self._mapping_id(mapping)]) /
                   float(dict(mapping.tags)["NH"])) # no. of mappings of read

    def _fraction_norm_by_mapping(self, mapping):
        return(1.0 / float(dict(mapping.tags)["NH"])) # no. of mappings of read

    def _fraction_norm_by_overlap(self, mapping):
        return(1.0 /
               float(self.mappings_and_no_of_overlaps[
                   self._mapping_id(mapping)]))

    def _overlapping_mappings(self, sam, entry):
        for mapping in sam.fetch(
                reference=entry.seq_id, start=entry.start, end=entry.end):
            if mapping.overlap(entry.start, entry.end) < self._min_overlap:
                continue
            yield(mapping)

    def _mapping_id(self, mapping):
        return(":".join([str(mapping.tid), mapping.qname, str(mapping.flag),
                         str(mapping.pos), str(mapping.aend)]))

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
        for annotation_file_path in sorted(path_and_name_combos.keys()):
            table_columns = []
            entries = []
            for entry in gff3_parser.entries(open(annotation_file_path)):
                entries.append(direction + "\t" + str(entry))
            table_columns.append(entries)
            for read_file, gene_quanti_path in path_and_name_combos[
                    annotation_file_path]:
                reader = csv.reader(open(gene_quanti_path), delimiter="\t")
                next(reader) # skip first line
                table_columns.append([row[column] for row in reader])
            # Generate a table by rotating the column list
            table = zip(*table_columns)
            for row in table:
                output_fh.write("\t".join(row) + "\n")
