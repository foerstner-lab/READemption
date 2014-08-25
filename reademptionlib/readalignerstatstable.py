class ReadAlignerStatsTable(object):

    def __init__(self, read_processing_stats, alignment_stats, 
                 primary_read_aligner_stats, realigner_stats, libs, output_path,
                 paired_end):
        self._table = []
        self._read_processing_stats = read_processing_stats
        self._alignment_stats = alignment_stats
        self._primary_read_aligner_stats = primary_read_aligner_stats
        self._realigner_stats = realigner_stats
        self._libs = libs
        self._output_path = output_path
        self._paired_end = paired_end
    
    def write(self):
        self._add_global_countings()
        self._add_reference_wise_coutings()
        with open(self._output_path, "w") as table_fh:
            table_fh.write("\n".join(["\t".join([str(cell) for cell in row])
                                      for row in self._table]) + "\n")

    def _add_global_countings(self):
        for title, data in [
            ("Libraries", self._libs),
            ("No. of input reads", 
             self._get_read_process_numbers("total_no_of_reads")),
            ("No. of reads - PolyA detected and removed", 
             self._get_read_process_numbers("polya_removed")),
            ("No. of reads - Single 3' A removed", 
             self._get_read_process_numbers("single_a_removed")),
            ("No. of reads - Unmodified", 
             self._get_read_process_numbers("unmodified")),
            ("No. of reads - Removed as too short", 
             self._get_read_process_numbers("too_short")),
            ("No. of reads - Long enough and used for alignment", 
             self._get_read_process_numbers("long_enough")),
            ("Total no. of aligned reads", 
             self._total_alignment_stat_numbers("no_of_aligned_reads")),
            ("Total no. of unaligned reads", 
             self._total_alignment_stat_numbers("no_of_unaligned_reads")),
            ("Total no. of uniquely aligned reads", 
             self._total_alignment_stat_numbers(
                    "no_of_uniquely_aligned_reads")),
            ("Total no. of alignments", 
             self._total_alignment_stat_numbers("no_of_alignments")),
            ("Total no. of split alignments", 
             self._total_alignment_stat_numbers("no_of_split_alignments")),
            ("Percentage of aligned reads (compared to no. of input reads)",
             self._perc_aligned_reads_all_input()),
            ("Percentage of aligned reads (compared to no. of long enough reads)",
             self._perc_aligned_reads_all_long_enough()),
            ("Percentage of uniquely aligned reads (in relation to all aligned "
             "reads)", self._perc_uniquely_aligned_reads())]:
            self._table.append([title] + data)

    def _add_reference_wise_coutings(self):
        ref_ids = sorted(list(list(self._alignment_stats.values())[0][
                    "stats_per_reference"].keys()))
        for ref_id in ref_ids:
            for title_template, data in [
                ("%s - No. of aligned reads", 
                 self._alignment_number_per_ref_seq(
                        ref_id, "no_of_aligned_reads")),
                ("%s - No. of uniquely aligned reads", 
                 self._alignment_number_per_ref_seq(
                        ref_id, "no_of_uniquely_aligned_reads")),
                ("%s - No. of alignments", 
                 self._alignment_number_per_ref_seq(
                        ref_id, "no_of_alignments")),
                ("%s - No. of split alignments", 
                 self._alignment_number_per_ref_seq(
                        ref_id, "no_of_split_alignments"))]:
                    self._table.append([title_template % ref_id] + data)

    def _alignment_number_per_ref_seq(self, ref_id, attribute):
        return [round(self._alignment_stats[lib]["stats_per_reference"][
                ref_id].get(attribute, 0)) for lib in self._libs]

    def _total_alignment_stat_numbers(self, attribute, round_nums=True):
        countings = [self._alignment_stats[lib]["stats_total"].get(
                attribute, 0) for lib in self._libs]
        if round_nums is True:
            return [round(counting) for counting in countings]
        else:
            return countings

    def _get_read_process_numbers(self, attribute):
        factor = 1
        if self._paired_end:
            factor = 2
        return [self._read_processing_stats[lib][attribute] * factor
                for lib in self._libs]

    def _perc_aligned_reads_all_input(self):
        return [
            round(self._calc_percentage(aligned_reads, total_reads), 2)
            for aligned_reads, total_reads in
            zip(self._total_alignment_stat_numbers(
                    "no_of_aligned_reads", round_nums=False),
                self._get_read_process_numbers("total_no_of_reads"))]

    def _perc_aligned_reads_all_long_enough(self):
        return [
            round(self._calc_percentage(aligned_reads, total_reads), 2)
            for aligned_reads, total_reads in
            zip(self._total_alignment_stat_numbers(
                    "no_of_aligned_reads", round_nums=False),
                self._get_read_process_numbers("long_enough"))]
    
    def _perc_uniquely_aligned_reads(self):
        return [
            round(self._calc_percentage(uniquely_aligned_reads, aligned_reads)
                  , 2)
            for uniquely_aligned_reads, aligned_reads  in zip(
                self._total_alignment_stat_numbers(
                    "no_of_uniquely_aligned_reads"),
                self._total_alignment_stat_numbers(
                    "no_of_aligned_reads"))]

    def _calc_percentage(self, mult, div):
        try:
            return float(mult)/float(div)*100
        except ZeroDivisionError:
            return 0.0
