import os
import concurrent.futures
from reademptionlib.paths import Paths
from reademptionlib.wiggle import WiggleWriter


class Controller_TK(object):

    def __init__(self, args):
        """Create an instance."""
        self._args = args
        self._paths = Paths

    def viz_align_TK(self):
        """Generate plots based on the read processing and mapping"""
        from reademptionlib.vizalign2 import AlignViz2
        align_viz_tk = AlignViz2()
        if self._args.input_align:
            align_viz_tk.alignment_stats(
                str(self._args.input_align),
                str(self._args.output_align),
                str(self._args.output_align))
        if self._args.input_process:
            align_viz_tk.process_stats(
                str(self._args.input_process),
                str(self._args.output_process),
                str(self._args.output_process))
        
    def get_and_submit_bam_stats(self):
        from reademptionlib.readalignerstats import ReadAlignerStats
        get_bam_stats = ReadAlignerStats()
        total_aligned_reads = []
        total_uniquely_aligned_reads = []
        all_aligned_reads = []
        for sample_or_stats_total, sub_dic in (
                get_bam_stats.count(self._args.input_coverage, "NA")).items():
            for stats_total, value in sub_dic.items():
                if sample_or_stats_total == 'stats_total':
                    if stats_total == 'no_of_uniquely_aligned_reads':
                        total_uniquely_aligned_reads.append(value)
                    if stats_total == 'no_of_aligned_reads':
                        total_aligned_reads.append(value)
                else:
                    for align_type, Nr in value.items():
                        if align_type == 'no_of_aligned_reads':
                            all_aligned_reads.append(Nr)
        min_no_of_aligned_reads = min(all_aligned_reads)
        if not self._args.normalize_by_uniquely:
            no_of_aligned_reads = "".join([
                '{:.1f}'.format(x) for x in total_aligned_reads])
        else:
            no_of_aligned_reads = "".join([
                '{:.1f}'.format(x) for x in total_uniquely_aligned_reads])
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=self._args.processes) as executor:
            jobs.append(executor.submit(
                self.create_coverage_files_TK, float(no_of_aligned_reads),
                float(min_no_of_aligned_reads)))
        self._check_job_completeness(jobs)

    def create_coverage_files_TK(
            self, no_of_aligned_reads, min_no_of_aligned_reads):
        if not self._args.non_strand_specific:
            strands = ["forward", "reverse"]
        else:
            strands = ["forward_and_reverse"]
        if self._args.skip_read_count_splitting:
            read_count_splitting = False
        else:
            read_count_splitting = True
        from reademptionlib.coveragecalculator import CoverageCalculator
        coverage_calculator = CoverageCalculator(
            read_count_splitting=read_count_splitting,
            uniquely_aligned_only=self._args.unique_only,
            coverage_style=self._args.coverage_style,
            clip_length=self._args.clip_length,
            non_strand_specific=self._args.non_strand_specific)
        (coverage_writers_raw, coverage_writers_tnoar_min_norm,
         coverage_writers_tnoar_mil_norm) = self._wiggle_writers(
             strands, no_of_aligned_reads, min_no_of_aligned_reads)
        for ref_seq, coverages in coverage_calculator.ref_seq_and_coverages(
                self._args.input_coverage):
            for strand in strands:
                coverage_writers_raw[strand].write_replicons_coverages(
                    ref_seq, coverages[strand])
            for strand in strands:
                coverage_writers_tnoar_min_norm[
                    strand].write_replicons_coverages(
                    ref_seq, coverages[strand],
                    factor=min_no_of_aligned_reads/no_of_aligned_reads)
            for strand in strands:
                coverage_writers_tnoar_mil_norm[
                    strand].write_replicons_coverages(
                    ref_seq, coverages[strand],
                    factor=self._args.factor/no_of_aligned_reads)
        for strand in strands:
            coverage_writers_raw[strand].close_file()

    def _wiggle_writers(self, strands, no_of_aligned_reads,
                        min_no_of_aligned_reads):
        if self._args.output_coverage is None:
            filename = self.get_filename(self._args.input_coverage)
            if os.path.dirname(self._args.input_coverage) == "":
                output_path = "."
            else:
                output_path = os.path.dirname(self._args.input_coverage)
        else:
            filename = self.get_filename(self._args.output_coverage)
            if os.path.dirname(self._args.output_coverage) == "":
                output_path = "."
            else:
                output_path = os.path.dirname(self._args.output_coverage)
        if self._args.trackname is None:
            trackname = filename
        else:
            trackname = self.get_filename(self._args.trackname)
        coverage_writers_raw = dict([(
            strand, WiggleWriter(
                "%s_%s" % (trackname, strand),
                open("%s/%s_%s-raw.wig" % (output_path, filename, strand),
                     "w"))) for strand in strands])
        coverage_writers_tnoar_min_norm = dict([(
            strand, WiggleWriter(
                "%s_%s" % (trackname, strand),
                open(
                    "%s/%s_%s-tnoar_mil_norm_multi_by_%.1f_div_by_%.1f.wig" % (
                        output_path, filename, strand,
                        min_no_of_aligned_reads, no_of_aligned_reads),
                    "w"))) for strand in strands])
        coverage_writers_tnoar_mil_norm = dict([(
            strand, WiggleWriter(
                "%s_%s" % (trackname, strand),
                open(
                    "%s/%s_%s-tnoar_mil_norm_multi_by_%.1f_div_by_%.1f.wig" % (
                        output_path, filename, strand, self._args.factor,
                        no_of_aligned_reads),
                    "w"))) for strand in strands])
        return(coverage_writers_raw, coverage_writers_tnoar_min_norm,
               coverage_writers_tnoar_mil_norm)

    def get_filename(self, path_file):
        return(os.path.splitext(os.path.basename(path_file))[0])

    def _check_job_completeness(self, jobs):
        """Check the completness of each job in a list"""
        for job in concurrent.futures.as_completed(jobs):
            if job.exception():
                raise(job.exception())
            


