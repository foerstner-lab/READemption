import sys
sys.path.append(".")
from libs.sam import SamParser

class GRCreator(object):

    def create_gr_files(self, read_file_names, read_mapping_result_paths, 
                        ref_ids_to_file_name, gr_folder):
        for read_file_name, read_mapping_result_path in zip(
            read_file_names, read_mapping_result_paths):
             for ref_seq_id, ref_seq_file_name in ref_ids_to_file_name.items():
                plus_strand_output_file = self._output_file_path(
                    gr_folder, read_file_name, ref_seq_file_name, "+")
                minus_strand_output_file = self._output_file_path(
                    gr_folder, read_file_name, ref_seq_file_name, "-")
                gr_file_builder = GRFileBuilder(
                    read_mapping_result_path, ref_seq_id, 
                    plus_strand_output_file, minus_strand_output_file)
                gr_file_builder.build_gr_files()

    def create_read_normalized_gr_files(
        self, read_file_names, read_mapping_result_paths, 
        ref_ids_to_file_name, gr_read_normalized_folder):
        read_file_names_and_mapped_reads = {}
        sam_parser = SamParser()
        for read_file_name, read_mapping_result_path in zip(
            read_file_names, read_mapping_result_paths):
            ref_seqs_and_mappings, ref_seqs_and_mapped_reads = (
                sam_parser.mapping_countings(open(read_mapping_result_path)))
            # Sum the number of mapped reads for all reference sequences
            read_file_names_and_mapped_reads[read_file_name] = sum(
                ref_seqs_and_mapped_reads.values())
        min_no_of_reads = min(read_file_names_and_mapped_reads.values())

        for read_file_name, read_mapping_result_path in zip(
            read_file_names, read_mapping_result_paths):
            norm_value = read_file_names_and_mapped_reads[read_file_name]
            for ref_seq_id, ref_seq_file_name in ref_ids_to_file_name.items():
                plus_strand_output_file = self._read_normalized_file_path(
                    gr_read_normalized_folder, read_file_name, 
                    ref_seq_file_name, "+", norm_value, min_no_of_reads)
                minus_strand_output_file = self._read_normalized_file_path(
                    gr_read_normalized_folder, read_file_name, 
                    ref_seq_file_name, "-", norm_value, min_no_of_reads)
                gr_file_builder = GRFileBuilder(
                    read_mapping_result_path, ref_seq_id, 
                    plus_strand_output_file, minus_strand_output_file,
                    normalization_value=norm_value, 
                    multiplier=min_no_of_reads)
                gr_file_builder.build_gr_files()    
        
    def _output_file_path(self, folder_path, 
                          read_file_name, ref_seq_file_name, strand):
        strand_string = {"-" : "minus", "+" : "plus"}[strand]
        return("%s/%s_in_%s.%s_strand.gr" % (
                folder_path, read_file_name, ref_seq_file_name, strand_string))

    def _read_normalized_file_path(
        self, folder_path, read_file_name, ref_seq_file_name, strand, 
        normalization_value, multiplier):
        strand_string = {"-" : "minus", "+" : "plus"}[strand]
        normalization_value = round(normalization_value, 2)
        multiplier = round(multiplier, 2)
        return("%s/%s_in_%s_norm_by_%s_mult_by_%s.%s_strand.gr" % (
                folder_path, read_file_name, ref_seq_file_name, 
                normalization_value, multiplier, strand_string))

class GRFileBuilder(object):

    def __init__(self, input_sam_path, ref_seq_id, plus_strand_output_file, 
                 minus_strand_output_file, normalization_value=1, 
                 multiplier=1):
        self.input_sam_path = input_sam_path
        self.ref_seq_id = ref_seq_id
        self.plus_strand_output_file = plus_strand_output_file
        self.minus_strand_output_file = minus_strand_output_file
        self.normalization_value = normalization_value
        self.multiplier = multiplier

    def build_gr_files(self):
        coverage_plus_strand, coverage_minus_strand = self._calc_raw_coverages(
            open(self.input_sam_path))
        if self._norm_or_multi_needed:
            coverage_plus_strand = self._normalize_and_multiply(
                coverage_plus_strand)
            coverage_minus_strand = self._normalize_and_multiply(
                coverage_minus_strand)
        self._build_gr_file(coverage_plus_strand, open(
                self.plus_strand_output_file, "w"))
        self._build_gr_file(coverage_minus_strand, open(
                self.minus_strand_output_file, "w"))
        
    def _norm_or_multi_needed(self):
        return(not(self.normalization_value == 1 and self.multiplier == 1))

    def _normalize_and_multiply(self, coverages):
        return([coverage * self.multiplier / self.normalization_value 
                for coverage in coverages])

    def _calc_raw_coverages(self, sam_fh):
        coverages_plus_strand = []
        coverages_minus_strand = []
        for entry in self._sam_entries(sam_fh):
            if not entry.reference == self.ref_seq_id:
                continue
            if entry.strand == "+":
                self._add_coverage(entry, coverages_plus_strand)
            elif entry.strand == "-":
                self._add_coverage(entry, coverages_minus_strand)
        coverages_minus_strand = [
            -1.0 * coverage for coverage in coverages_minus_strand]
        return(coverages_plus_strand, coverages_minus_strand)

    def _sam_entries(self, sam_fh):
        sam_parser = SamParser()
        for entry in sam_parser.entries(sam_fh):
            yield(entry)

    def _add_coverage(self, entry, coverages):
        if len(coverages) < entry.end:
            self._extend_coverages(coverages, entry.end)
        for position in range(entry.start-1, entry.end):
            coverages[position] += 1.0 / float(entry.number_of_hits_as_int)  
            
    def _extend_coverages(self, coverages, end):
        coverages.extend([0] * (end - len(coverages)))

    def _build_gr_file(self, coverages, output_fh):
        for pos, coverage in enumerate(coverages):
            # Skip zero values to save file space
            if coverage == 0.0:
                continue
            output_fh.write("%s\t%s\n" % (pos, coverage))

