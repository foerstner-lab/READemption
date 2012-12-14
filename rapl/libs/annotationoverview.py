from libs.intersectbed import IntersectBedParser
from libs.gff3 import Gff3Parser
from libs.sam import SamParser

class AnnotationOverview(object):

    def __init__(self, min_overlap=10):
        self._annotation_hit_countings = {}
        self._mappings_per_read = {}
        self._mappings_and_no_of_overlaps = {}
        self._intersect_bed_parser = IntersectBedParser()
        self._gff3_parser = Gff3Parser()
        self._sam_parser = SamParser()
        self._sense_str = "s"
        self._antisense_str = "a"
        self._min_overlap = min_overlap

    def init_counting_table(
        self, read_file_name, annotation_file, annotation_file_path):
        """Initiate the counting table with 0 values

        This table will hold the values for all read mapping /
        annotation combination.
        """
        self._annotation_hit_countings.setdefault(read_file_name, {})
        self._annotation_hit_countings[read_file_name].setdefault(
            annotation_file, {})
        for entry in self._gff3_parser.entries(open(annotation_file_path)):
            key = self._gff_entry_to_gene_key(entry)
            self._annotation_hit_countings[read_file_name][annotation_file][key] = {
                self._sense_str : 0, self._antisense_str :  0}

    def _gff_entry_to_gene_key(self, entry):
        """Generate a dictionary key based on an GFF entry.

        The keys should have the same style as the ones produced by
        _intersectbed_entry_to_gene_key.

        """
        return(self._values_to_gene_key(entry.seq_id, entry.feature, entry.start,
                                   entry.end, entry.strand))

    def _intersectbed_entry_to_gene_key(self, entry):
        """Generate a dictionary key based on an intersectBed entry.

        The keys should have the same style as the ones produced by
        _gff_entry_to_gene_key

        """
        return(self._values_to_gene_key(
                entry.gff_seq_id, entry.gff_feature, entry.gff_start,
                entry.gff_end, entry.gff_strand))

    def _values_to_gene_key(self, seq_id, feature, start, end, strand):
        return("|".join(
                [str(val) for val in [seq_id, feature, start, end, strand]]))

    def _intersectbed_entry_to_mapping_key(self, entry):
        return("|".join([
                    str(val) for val in[
                        entry.sam_reference, entry.sam_strand, entry.sam_end,
                        entry.sam_query_id]]))

    def get_read_mapping_freq(self, bam_file_path):
        for entry in self._sam_parser.entries_bam(bam_file_path):
            self._mappings_per_read[entry.query_id] = int(
                entry.number_of_hits_as_int)

        # get_overlapping_genes_of_mappings
    def get_mapping_overlap_freq(
        self, read_file_name, annotation_hit_file_path):
        """Store the genes each mapping is overlapping with

        This will be later used to count the number of genes each
        mapping is overlapping with.

        """
        self._mappings_and_no_of_overlaps[read_file_name] = {}
        for entry in self._intersect_bed_parser.entries(
            open(annotation_hit_file_path)):
            if not self._valid_overlap(entry):
                continue
            mapping_key = self._intersectbed_entry_to_mapping_key(entry)
            gene_key = self._intersectbed_entry_to_gene_key(entry)
            self._mappings_and_no_of_overlaps[read_file_name].setdefault(
                mapping_key, set())
            self._mappings_and_no_of_overlaps[
                read_file_name][mapping_key].add(gene_key)

    def _valid_overlap(self, entry):
        return(entry.overlap >= self._min_overlap)

    def count_overlapping_reads_per_gene(
        self, read_file_name, annotation_file, annotation_hit_file_path):
        for entry in self._intersect_bed_parser.entries(
            open(annotation_hit_file_path)):
            if not self._valid_overlap(entry):
                continue
            gene_key = self._intersectbed_entry_to_gene_key(entry)
            mapping_key = self._intersectbed_entry_to_mapping_key(entry)
            fraction = (
                float(1)
                / float(self._mappings_per_read[entry.sam_query_id])
                / float(len(self._mappings_and_no_of_overlaps[
                    read_file_name][mapping_key])))
            self._annotation_hit_countings[
                read_file_name][annotation_file][gene_key][
                entry.orientation] += fraction

    def clean_dicts(self):
        self._mappings_per_read = {}
        self._mappings_and_no_of_overlaps = {}

    def write_overview_tables(
        self, annotation_file, annotation_file_path, read_file_names,
        output_file_sense, output_file_antisense,
        total_numbers_of_mapped_reads=None):
        output_sense_fh = open(output_file_sense, "w")
        output_antisense_fh = open(output_file_antisense, "w")
        self._write_header(output_sense_fh, read_file_names,
                           total_numbers_of_mapped_reads)
        self._write_header(output_antisense_fh, read_file_names)
        for entry in self._gff3_parser.entries(open(annotation_file_path)):
            gene_key = self._gff_entry_to_gene_key(entry)
            # Write sense output line
            output_sense_fh.write(str(entry) + "\t")
            output_sense_fh.write(
                "\t".join([str(self._annotation_hit_countings[read_file_name][
                            annotation_file][gene_key][self._sense_str])
                     for read_file_name in read_file_names]) + "\n")
            # Write anti-sense output line
            output_antisense_fh.write(str(entry) + "\t")
            output_antisense_fh.write(
                "\t".join([str(self._annotation_hit_countings[read_file_name][
                            annotation_file][gene_key][self._antisense_str])
                     for read_file_name in read_file_names]) + "\n")

    def _write_header(
        self, output_fh, read_file_names, total_numbers_of_mapped_reads=None):
        if total_numbers_of_mapped_reads:
            print(total_numbers_of_mapped_reads)
            output_fh.write(
                "\t".join([""] * 9 +
                          ["%s (normalized by %s mapped reads)" % (
                            read_file_name, total_number_of_mapped_reads)
                           for read_file_name, total_number_of_mapped_reads
                           in  zip(read_file_names, total_numbers_of_mapped_reads)])
                + "\n")
        else:
            output_fh.write("\t".join([""] * 9 + read_file_names) + "\n")


