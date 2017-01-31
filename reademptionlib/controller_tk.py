import os
import concurrent.futures
from reademptionlib.helpers import Helpers
from reademptionlib.paths import Paths
from reademptionlib.readalignerstats import ReadAlignerStats
from reademptionlib.rawstatdata import RawStatDataWriter
from reademptionlib.wiggle import WiggleWriter


class Controller_TK(object):

    def __init__(self, args):
        """Create an instance."""
        self._args = args
        self._paths = Paths(args)
        self._helpers = Helpers(args)

    def generate_read_alignment_stats(
            self, filenames, BAM_files, output_path, output_file):
        """Manage the generation of alingment statistics."""
        raw_stat_data_writer = RawStatDataWriter(pretty=True)
        read_files_and_jobs = {}
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=self._args.processes) as executor:
            for (filename, BAM_file) in zip(
                    filenames, BAM_files):
                read_aligner_stats = ReadAlignerStats()
                read_files_and_jobs[filename] = executor.submit(
                    read_aligner_stats.count, BAM_files, "NA")
        # Evaluate thread outcome
        self._helpers.check_job_completeness(read_files_and_jobs.values())
        read_files_and_stats = dict(
            [(filename, job.result())
             for filename, job in read_files_and_jobs.items()])
        raw_stat_data_writer.write(
            read_files_and_stats, output_file)

    def viz_align_tk(self):
        """Generate plots based on the read processing and mapping"""
        from reademptionlib.vizalign import AlignViz
        align_viz_tk = AlignViz()
        if self._args.input_BAM:
            for (input_file, output_file) in zip(
                    self._args.input_BAM, self._args.output_BAM):
                self.generate_read_alignment_stats(
                    self.get_filename(
                        input_file), input_file,
                    os.path.dirname(output_file), output_file)
        # Missing: the automated mapping process
        if self._args.input_process and not self._args.input_align:
            align_viz_tk.processing_viz(self._args.input_process,
                                        self._args.output_path)
        if self._args.input_align and not self._args.input_process:
            align_viz_tk.alignment_viz(self._args.input_align,
                                       self._args.output_path)
        if self._args.input_process and self._args.input_align:
            align_viz_tk.processing_viz(self._args.input_process,
                                        self._args.output_path)
            align_viz_tk.alignment_viz(self._args.input_align,
                                       self._args.output_path)
            align_viz_tk.alignment_processing_overview(
                self._args.input_process, self._args.input_align,
                self._args.output_path)

    def viz_deseq_tk(self):
        from reademptionlib.vizdeseq import DESeqViz
        deseq_viz_tk = DESeqViz(self._args.input_file, self._args.output_path,
                                self._args.padj_cutoff, self._args.condition,
                                self._args.alpha, self._args.color_sig,
                                self._args.color_non_sig, self._args.shape,
                                self._args.glyph_size)
        deseq_viz_tk.read_and_modificate_input()

    def viz_gene_quanti_tk(self):
        from reademptionlib.vizgenequanti import GeneQuantiViz
        gene_quanti_viz_tk = GeneQuantiViz(self._args.input_file,
                                           self._args.lib_names,
                                           self._args.output_path)
        gene_quanti_viz_tk.parse_input_table()

    def get_and_submit_bam_stats(self):
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
        self._helpers.check_job_completeness(jobs)

    def create_coverage_files_tk(
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


