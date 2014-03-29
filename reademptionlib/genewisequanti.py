import csv
import os.path
from reademptionlib.gff3 import Gff3Parser
import pysam

class GeneWiseQuantification(object):

    def __init__(self, min_overlap=1, norm_by_alignment_freq=True,
                 norm_by_overlap_freq=True, allowed_features_str=None,
                 skip_antisense=False, unique_only=False):
        """
        - normalize_by_alignment: consider that some reads are aligned at
          more than one location and only count fractions
        - normalize_by_overlapping_genes: consider that some alignment
          overlap with more than on gene

        """
        self._min_overlap = min_overlap
        self._norm_by_alignment_freq = norm_by_alignment_freq
        self._norm_by_overlap_freq = norm_by_overlap_freq
        self._allowed_features = _allowed_features(allowed_features_str)
        self._skip_antisense = skip_antisense
        self._unique_only = unique_only

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
                if _entry_to_use(entry, self._allowed_features) is False:
                    continue
                for alignment in self._overlapping_alignments(sam, entry):
                    alignment_id = self._alignment_id(alignment)
                    self.alignments_and_no_of_overlaps.setdefault(alignment_id, 0)
                    self.alignments_and_no_of_overlaps[alignment_id] += 1

    def quantify(self, read_alignment_path, annotation_path, output_path,
                 pseudocounts=False):
        self._quantify(read_alignment_path, annotation_path, output_path,
                      self._fraction_calc_method(), pseudocounts)

    def _quantify(self, read_alignment_path, annotation_path, output_path,
                  fraction_calc_method, pseudocounts=False):
        sam = pysam.Samfile(read_alignment_path)
        gff3_parser = Gff3Parser()
        output_fh = open(output_path, "w")
        output_fh.write("#" + "\t".join(_gff_field_descriptions() 
                                        + ["sense", "antisense"]) + "\n")
        for entry in gff3_parser.entries(open(annotation_path)):
            if _entry_to_use(entry, self._allowed_features) is False:
                continue
            if pseudocounts is False:
                sum_sense = 0
                sum_antisense = 0
            else:
                sum_sense = 1
                sum_antisense = 1
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
        if alignment.is_read2 == False:
            if ((entry.strand == "+" and alignment.is_reverse is False) or
                (entry.strand == "-" and alignment.is_reverse is True)):
                return True
        # Mate pair for paired end sequencing
        elif alignment.is_read2 == True:
            if ((entry.strand == "+" and alignment.is_reverse is True) or
                (entry.strand == "-" and alignment.is_reverse is False)):
                return True
        return False 

    def _fraction_calc_method(self):
        if self._norm_by_alignment_freq and self._norm_by_overlap_freq:
           return self._fraction_norm_by_alignment_and_overlap
        elif self._norm_by_alignment_freq and not self._norm_by_overlap_freq:
            return self._fraction_norm_by_alignment
        elif not self._norm_by_alignment_freq and self._norm_by_overlap_freq:
            return self._fraction_norm_by_overlap
        return self._fraction_calc_constant_one

    def _alignment_tags(self, alignment):
        return dict(alignment.tags)

    def _fraction_calc_constant_one(self, alignment):
        return 1.0

    def _fraction_norm_by_alignment_and_overlap(self, alignment):
        alignment_tags = self._alignment_tags(alignment)
        return (1.0 /
                float(self.alignments_and_no_of_overlaps[
                    self._alignment_id(alignment)]) /
                float(alignment_tags["NH"]) / # no. of alignments of read
                float(alignment_tags.get("XL", 1))) # no. of splits

    def _fraction_norm_by_alignment(self, alignment):
        alignment_tags = self._alignment_tags(alignment)
        return (1.0 / float(alignment_tags["NH"]) / # no. of alignments of read
                float(alignment_tags.get("XL", 1))) # no. of splits

    def _fraction_norm_by_overlap(self, alignment):
        alignment_tags = self._alignment_tags(alignment)
        return (1.0 /
                float(self.alignments_and_no_of_overlaps[
                    self._alignment_id(alignment)]) /
                float(alignment_tags.get("XL", 1))) # no. of splits

    def _overlapping_alignments(self, sam, entry):
        # The substraction of 1 from the start is necessary to perform
        # this correctly (checked in IGB, IGV and the unit testings).
        for alignment in sam.fetch(
                reference=entry.seq_id, start=entry.start-1, end=entry.end):
            if alignment.overlap(entry.start-1, entry.end) < self._min_overlap:
                continue
            if self._skip_antisense:
                if not self._same_strand(entry, alignment):
                    continue
            if self._unique_only:
                if dict(alignment.tags)["NH"] != 1:
                    continue
            yield(alignment)

    def _alignment_id(self, alignment):
        return (":".join([str(alignment.tid), alignment.qname, str(alignment.flag),
                          str(alignment.pos), str(alignment.aend)]))

    def _values_to_gene_key(self, seq_id, feature, start, end, strand):
        return ("|".join(
                [str(val) for val in [seq_id, feature, start, end, strand]]))

class GeneWiseOverview(object):

    def __init__(self, allowed_features_str=None, skip_antisense=False):
        self._allowed_features = _allowed_features(allowed_features_str)
        self._skip_antisense = skip_antisense

    def create_overview_raw_countings(
            self, path_and_name_combos, read_files, overview_path):
        self._create_overview(path_and_name_combos, read_files, overview_path)

    def create_overview_rpkm(
            self, path_and_name_combos, read_files, overview_path, libs_and_tnoar):
        self._create_overview(path_and_name_combos, read_files, overview_path,
                              normalization="RPKM", libs_and_tnoar=libs_and_tnoar)

    def create_overview_norm_by_tnoar(
            self, path_and_name_combos, read_files, overview_path, libs_and_tnoar):
        self._create_overview(path_and_name_combos, read_files, overview_path,
                              normalization="TNOAR", libs_and_tnoar=libs_and_tnoar)

    def _create_overview(self, path_and_name_combos, read_files, overview_path,
                         normalization=None, libs_and_tnoar=None):
        output_fh = open(overview_path, "w")
        # Write header
        output_fh.write("\t".join(
                ["Orientation of counted reads relative to the strand "
                 "location of the annotation"] + _gff_field_descriptions() 
                + read_files) + "\n")
        self._add_to_overview(
            path_and_name_combos, "sense", 9, output_fh, normalization,
            libs_and_tnoar)
        if self._skip_antisense is False:
            self._add_to_overview(
                path_and_name_combos, "anti-sense", 10, output_fh,
                normalization, libs_and_tnoar)

    def _add_to_overview(
            self, path_and_name_combos, direction, column, output_fh,
            normalization=None, libs_and_tnoar=None):
        gff3_parser = Gff3Parser()
        for annotation_path in sorted(path_and_name_combos.keys()):
            table_columns = []
            entries = []
            seq_lengths = []
            for entry in gff3_parser.entries(open(annotation_path)):
                if _entry_to_use(entry, self._allowed_features) is False:
                    continue
                entries.append(direction + "\t" + str(entry))
                seq_lengths.append(entry.end - entry.start + 1)
            table_columns.append(entries)
            for read_file, gene_quanti_path in path_and_name_combos[
                    annotation_path]:
                reader = csv.reader(open(gene_quanti_path), delimiter="\t")
                next(reader) # skip first line
                if normalization == "RPKM":
                    table_columns.append([
                        self._rpkm(row[column], length, libs_and_tnoar[read_file])
                        for row, length in zip(reader, seq_lengths)])
                elif normalization == "TNOAR":
                    table_columns.append([
                        self._norm_by_tnoar(row[column], libs_and_tnoar[read_file])
                        for row, length in zip(reader, seq_lengths)])
                else:
                    table_columns.append([row[column] for row in reader])
            # Generate a table by rotating the column list
            table = zip(*table_columns)
            for row in table:
                output_fh.write("\t".join(row) + "\n")

    def _rpkm(self, counting, length, total_no_of_aligned_reads):
        """
        Formula in Supplemenatary Material S1 of
        http://www.nature.com/nmeth/journal/v5/n7/full/nmeth.1226.html

        R = (10^9 * C) / (N * L)

        with C = is the number of mappable reads that fell onto the gene 
             N = total number of mappable read
             L = length of the gene
        
        """
        return str(float(counting)*float(10**9) /
                   (float(total_no_of_aligned_reads)*float(length)))

    def _norm_by_tnoar(self, counting, total_no_of_aligned_reads):
        return str(float(counting)/float(total_no_of_aligned_reads))

def _entry_to_use(entry, allowed_features):
    if allowed_features is None:
        return True
    if entry.feature in allowed_features:
        return True
    return False 

def _allowed_features(allowed_features_str):
    if allowed_features_str is None:
        return None
    else:
        return [
        feature.strip() for feature in allowed_features_str.split(",")]

def _gff_field_descriptions():
    return ["Sequence name", "Source", "Feature", "Start", "End", "Score", 
            "Strand", "Frame", "Attributes"]
