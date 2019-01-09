import concurrent.futures
import sys
from reademptionlib.helpers import Helpers
from reademptionlib.paths import Paths
from reademptionlib.coveragecalculator import CoverageCalculator
from reademptionlib.rawstatdata import RawStatDataReader
from reademptionlib.wiggle import WiggleWriter


class CoverageController(object):

    def __init__(self, args):
        self._args = args
        self._paths = Paths(args)
        self._helpers = Helpers(args)

    def create_coverage_files(self):
        """Create coverage files based on the read alignments.
        The coverages are calculated per replicon and the results are
        written to the output file. This might be slower but if all
        coverages are detmined at once the data structure will become
        too large when working with large reference sequences.
        """
        self._helpers.test_folder_existance(
            self._paths.required_coverage_folders())
        raw_stat_data_reader = RawStatDataReader()
        alignment_stats = [raw_stat_data_reader.read(
            self._paths.read_alignments_stats_path)]
        lib_names = list(alignment_stats[0].keys())
        was_paired_end_alignment = self._helpers.was_paired_end_alignment(
            lib_names)
        if not was_paired_end_alignment:
            self._paths.set_read_files_dep_file_lists_single_end(
                self._paths.get_read_files(), lib_names)
        else:
            self._paths.set_read_files_dep_file_lists_paired_end(
                self._paths.get_read_files(), lib_names)
        # Get number of aligned or number of uniquely aligned reads
        if not self._args.normalize_by_uniquely:
            aligned_counting = "no_of_aligned_reads"
        else:
            aligned_counting = "no_of_uniquely_aligned_reads"
        read_files_aligned_read_freq = dict([
            (read_file,
             round(attributes["stats_total"][aligned_counting]))
            for read_file, attributes in alignment_stats[0].items()])
        min_no_of_aligned_reads = float(min(
            read_files_aligned_read_freq.values()))
        # Run the generation of coverage in parallel
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=self._args.processes) as executor:
            for lib_name, bam_path in zip(
                    lib_names, self._paths.read_alignment_bam_paths):
                no_of_aligned_reads = float(
                    read_files_aligned_read_freq[lib_name])
                jobs.append(executor.submit(
                    self._create_coverage_files_for_lib,
                    lib_name, bam_path, no_of_aligned_reads,
                    min_no_of_aligned_reads))
        # Evaluate thread outcome
        self._helpers.check_job_completeness(jobs)

    def _create_coverage_files_for_lib(
        self, lib_name, bam_path, no_of_aligned_reads,
            min_no_of_aligned_reads):
        """Perform the coverage calculation for a given library."""
        if not self._args.non_strand_specific:
            strands = ["forward", "reverse"]
        else:
            strands = ["forward_and_reverse"]
        if self._all_coverage_file_exist(
                lib_name, strands, no_of_aligned_reads,
                min_no_of_aligned_reads):
            return
        read_count_splitting = True
        if self._args.skip_read_count_splitting:
            read_count_splitting = False
        coverage_calculator = CoverageCalculator(
            read_count_splitting=read_count_splitting,
            uniquely_aligned_only=self._args.unique_only,
            coverage_style=self._args.coverage_style,
            clip_length=self._args.clip_length,
            non_strand_specific=self._args.non_strand_specific)
        (coverage_writers_raw, coverage_writers_tnoar_min_norm,
         coverage_writers_tnoar_mil_norm) = self._wiggle_writers(
            lib_name, strands, no_of_aligned_reads, min_no_of_aligned_reads)
        for ref_seq, coverages in coverage_calculator.ref_seq_and_coverages(
                bam_path):
            for strand in strands:
                coverage_writers_raw[strand].write_replicons_coverages(
                    ref_seq, coverages[strand])
                coverage_writers_tnoar_min_norm[
                    strand].write_replicons_coverages(
                    ref_seq, coverages[strand],
                    factor=min_no_of_aligned_reads/no_of_aligned_reads)
                coverage_writers_tnoar_mil_norm[
                    strand].write_replicons_coverages(
                    ref_seq, coverages[strand],
                    factor=1000000/no_of_aligned_reads)
        for strand in strands:
            coverage_writers_raw[strand].close_file()

    def _all_coverage_file_exist(
        self, lib_name, strands, no_of_aligned_reads,
            min_no_of_aligned_reads):
        """Test the existance of all coverage file of a library"""
        files = []
        for strand in strands:
            files.append(self._paths.wiggle_file_raw_path(lib_name, strand))
            files.append(self._paths.wiggle_file_tnoar_norm_min_path(
                lib_name, strand, multi=min_no_of_aligned_reads,
                div=no_of_aligned_reads))
            files.append(self._paths.wiggle_file_tnoar_norm_mil_path(
                lib_name, strand, multi=1000000,
                div=no_of_aligned_reads))
        if not any([self._helpers.file_needs_to_be_created(file, quiet=True)
                    for file in files]):
            sys.stderr.write(
                "The files %s exists. Skipping their generation.\n" %
                ", " .join(files))
            return True
        return False

    def _wiggle_writers(self, lib_name, strands, no_of_aligned_reads,
                        min_no_of_aligned_reads):
        """Write the calculated coverages to wiggle files."""
        coverage_writers_raw = dict([(
            strand, WiggleWriter(
                "%s_%s" % (lib_name, strand),
                open(self._paths.wiggle_file_raw_path(lib_name, strand),
                     "w"))) for strand in strands])
        coverage_writers_tnoar_min_norm = dict([(
            strand, WiggleWriter(
                "%s_%s" % (lib_name, strand),
                open(self._paths.wiggle_file_tnoar_norm_min_path(
                    lib_name, strand, multi=min_no_of_aligned_reads,
                    div=no_of_aligned_reads), "w"))) for strand in strands])
        coverage_writers_tnoar_mil_norm = dict([(
            strand, WiggleWriter(
                "%s_%s" % (lib_name, strand),
                open(self._paths.wiggle_file_tnoar_norm_mil_path(
                    lib_name, strand, multi=1000000,
                    div=no_of_aligned_reads), "w"))) for strand in strands])
        return (coverage_writers_raw, coverage_writers_tnoar_min_norm,
                coverage_writers_tnoar_mil_norm)
