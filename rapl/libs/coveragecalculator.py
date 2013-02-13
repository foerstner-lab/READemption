import pysam

class CoverageCalculator(object):

    def __init__(self, read_count_splitting=True, uniqueley_aligned_only=False,
                 first_base_only=False):
        self._read_count_splitting = read_count_splitting
        self._uniqueley_aligned_only = uniqueley_aligned_only
        self._first_base_only = first_base_only

    def ref_seq_and_coverages(self, bam_path):
        bam = self._open_bam_file(bam_path)
        self._coverage_add_function = self._select_coverage_add_function()
        self._coverages = {}
        for ref_seq, length in zip(bam.references, bam.lengths):
            for strand in ["forward", "reverse"]:
                self._coverages[strand] = [0.0] * length
            self._calc_coverage(ref_seq, bam)
            yield(ref_seq, self._coverages)

    def _calc_coverage(self, ref_seq, bam):
        for entry in bam.fetch(ref_seq):
            number_of_hits = dict(entry.tags)["NH"]
            if self._uniqueley_aligned_only is True and number_of_hits != 1:
                continue
            # Note: No translation from SAMParsers coordinates to python
            # list coorindates is needed.
            start = entry.pos
            end = entry.aend
            # Normalize coverage increment by number of read alignments
            # per read
            if self._read_count_splitting is True:
                increment = 1.0 / float(number_of_hits)
            else:
                increment = 1.0
            self._coverage_add_function(entry, increment, start, end)

    def _select_coverage_add_function(self):
        if self._first_base_only is False:
            return(self._add_whole_alignment_coverage)
        else:
            return(self._add_first_base_coverage)

    def _open_bam_file(self, bam_file):
        return(pysam.Samfile(bam_file))

    def _close_bam_fh(self, bam_fh):
        bam_fh.close()

    def _add_whole_alignment_coverage(self, entry, increment, start, end):
        if entry.is_reverse is False:
            self._coverages["forward"][start:end] = [
                coverage + increment for coverage in
                self._coverages["forward"][start:end]]
        else:
            self._coverages["reverse"][start:end] = [
                coverage - increment for coverage in
                self._coverages["reverse"][start:end]]

    def _add_first_base_coverage(self, entry, increment, start, end):
        if entry.is_reverse is False:
            self.replicons_and_coverages[
                "forward"][ref_id][start] = self.replicons_and_coverages[
                    "forward"][ref_id][start] + increment
        else:
            self.replicons_and_coverages[
                "reverse"][ref_id][end-1] = self.replicons_and_coverages[
                    "reverse"][ref_id][end-1] - increment
    
