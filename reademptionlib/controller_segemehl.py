import concurrent.futures

from reademptionlib.helpers import Helpers
from reademptionlib.paths import Paths
from reademptionlib.segemehl import Segemehl
from reademptionlib.readalignerstatstable import ReadAlignerStatsTable
from reademptionlib.rawstatdata import RawStatDataWriter, RawStatDataReader
from reademptionlib.readalignerstats import ReadAlignerStats
from reademptionlib.sambamconverter import SamToBamConverter


class SegemehlController():

    def __init__(self, args):
        self._args = args
        self._paths = Paths(args)
        self._helpers = Helpers(args)
        
    def align_with_segemehl(self):
        self._ref_seq_files = self._paths.get_ref_seq_files()
        self._paths.set_ref_seq_paths(self._ref_seq_files)
        if not self._args.paired_end:
            self._read_files = self._paths.get_read_files()
            self._lib_names = self._paths.get_lib_names_single_end()
            self._paths.set_read_files_dep_file_lists_single_end(
                self._read_files, self._lib_names)
            self._paths.set_read_files_dep_file_lists_single_end(
                self._read_files, self._lib_names)
            if not self._args.realign:
                self._set_primary_aligner_paths_to_final_paths()
            self._align_single_end_reads()
        self._sam_to_bam(
            self._paths.primary_read_aligner_sam_paths,
            self._paths.primary_read_aligner_bam_prefix_paths,
            self._paths.primary_read_aligner_bam_paths)
        self._generate_read_alignment_stats(
            self._lib_names,
            self._paths.primary_read_aligner_bam_paths,
            self._paths.unaligned_reads_paths,
            self._paths.primary_read_aligner_stats_path)
        final_unaligned_reads_paths = self._paths.unaligned_reads_paths
        if self._args.realign:
            self._run_realigner_and_process_alignments()
            self._merge_bam_files()
            final_unaligned_reads_paths = (
                self._paths.realigned_unaligned_reads_paths)
        if self._args.crossalign_cleaning_str is not None:
            self._remove_crossaligned_reads()

        """
        TODO
        Does currently not work as the stats for the read processing
        are not generated.
        self._generate_read_alignment_stats(
            self._lib_names,
            self._paths.read_alignment_bam_paths,
            final_unaligned_reads_paths,
            self._paths.read_alignments_stats_path)
        self._write_alignment_stat_table()
        self._align_viz.alignment_viz(
            self._paths.read_alignments_stats_path, "{}".format(
                self._paths.viz_align_base_folder))
        self._align_viz.processing_viz(
            self._paths.read_processing_stats_path, "{}".format(
                self._paths.viz_align_base_folder))
        self._align_viz.alignment_processing_overview(
            self._paths.read_processing_stats_path,
            self._paths.read_alignments_stats_path, "{}".format(
                self._paths.viz_align_base_folder))
        """
            
    def _set_primary_aligner_paths_to_final_paths(self):
        # If no remapping is performed the paths of the final bam files
        # is the paths of the primary mapper
        self._paths.primary_read_aligner_bam_prefix_paths = (
            self._paths.read_alignment_bam_prefix_paths)
        self._paths.primary_read_aligner_bam_paths = (
            self._paths.read_alignment_bam_paths)
        self._paths.primary_read_aligner_stats_path = (
            self._paths.read_alignments_stats_path)
            
    def _align_single_end_reads(self):
        """Manage the actual alignment of single end reads."""
        read_aligner = Segemehl(
            self._args.segemehl_bin, self._args.progress)
        if self._helpers.file_needs_to_be_created(self._paths.index_path):
            read_aligner.build_index(
                self._paths.ref_seq_paths, self._paths.index_path)
        for read_path, output_path, nomatch_path, bam_path in zip(
            self._paths.processed_read_paths,
            self._paths.primary_read_aligner_sam_paths,
            self._paths.unaligned_reads_paths,
                self._paths.read_alignment_bam_paths):
            if not self._helpers.file_needs_to_be_created(output_path):
                continue
            elif not self._helpers.file_needs_to_be_created(bam_path):
                continue
            read_aligner.run_alignment(
                read_path, self._paths.index_path, self._paths.ref_seq_paths,
                output_path, nomatch_path, int(self._args.processes),
                int(self._args.hit_strategy),
                int(self._args.segemehl_accuracy),
                float(self._args.segemehl_evalue), self._args.split,
                paired_end=False)
            
    def _sam_to_bam(self, sam_paths, bam_prefixes_paths, bam_paths):
        """Manage the conversion of mapped read from SAM to BAM format."""
        sam_to_bam_converter = SamToBamConverter()
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=self._args.processes) as executor:
            for sam_path, bam_prefix_path, bam_path in zip(
                    sam_paths, bam_prefixes_paths, bam_paths):
                if not self._helpers.file_needs_to_be_created(bam_path):
                    continue
                jobs.append(executor.submit(
                    sam_to_bam_converter.sam_to_bam,
                    sam_path, bam_prefix_path))
        # Evaluate thread outcome
        self._helpers.check_job_completeness(jobs)

    def _generate_read_alignment_stats(
            self, lib_names, result_bam_paths, unaligned_reads_paths,
            output_stats_path):
        """Manage the generation of alingment statistics."""
        raw_stat_data_writer = RawStatDataWriter(pretty=True)
        read_files_and_jobs = {}
        if not self._helpers.file_needs_to_be_created(output_stats_path):
            return
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=self._args.processes) as executor:
            for (lib_name, read_alignment_bam_path,
                 unaligned_reads_path) in zip(
                    lib_names, result_bam_paths, unaligned_reads_paths):
                read_aligner_stats = ReadAlignerStats()
                read_files_and_jobs[lib_name] = executor.submit(
                    read_aligner_stats.count, read_alignment_bam_path,
                    unaligned_reads_path)
        # Evaluate thread outcome
        self._helpers.check_job_completeness(read_files_and_jobs.values())
        read_files_and_stats = dict(
            [(lib_name, job.result())
             for lib_name, job in read_files_and_jobs.items()])
        raw_stat_data_writer.write(read_files_and_stats, output_stats_path)
        
    def _write_alignment_stat_table(self):
        """Manage the creation of the mapping statistic output table."""
        raw_stat_data_reader = RawStatDataReader()
        read_processing_stats = raw_stat_data_reader.read(
            self._paths.read_processing_stats_path)
        final_alignment_stats = raw_stat_data_reader.read(
            self._paths.read_alignments_stats_path)
        realignment_stats = None
        primary_aligner_stats = None
        if self._args.realign:
            primary_aligner_stats = raw_stat_data_reader.read(
                self._paths.primary_read_aligner_stats_path)
            realignment_stats = raw_stat_data_reader.read(
                self._paths.read_realigner_stats_path)
        read_aligner_stats_table = ReadAlignerStatsTable(
            read_processing_stats, final_alignment_stats,
            primary_aligner_stats,
            realignment_stats, self._lib_names,
            self._paths.read_alignment_stats_table_path, self._args.paired_end)
        read_aligner_stats_table.write()
