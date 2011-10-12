import sys
sys.path.append(".")
from libs.sam import SamParser

class GRCreator(object):
    pass

#ref_ids_to_file_name, 

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
        coverage_plus_strand, coverage_plus_minus = self._calc_raw_coverages(
            open(input_sam_path))
        if self._norm_or_multi_needed:
            coverage_plus_strand = self._normalize_and_multiply(
                coverage_plus_strand)
            coverage_minus_strand = self._normalize_and_multiply(
                coverage_minus_strand)
            
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

    def _build_gr_file(self, ref_seq_id, strand, output_fh):
        pass
    
