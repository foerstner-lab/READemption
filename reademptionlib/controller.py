import concurrent.futures
import os
import sys
import json
from reademptionlib.bammerger import BamMerger
from reademptionlib.coveragecalculator import CoverageCalculator
from reademptionlib.crossalignfilter import CrossAlignFilter
from reademptionlib.deseq import DESeqRunner
from reademptionlib.fasta import FastaParser
from reademptionlib.genewisequanti import GeneWiseQuantification
from reademptionlib.genewisequanti import GeneWiseOverview
from reademptionlib.paths import Paths
from reademptionlib.projectcreator import ProjectCreator
from reademptionlib.rawstatdata import RawStatDataWriter, RawStatDataReader
from reademptionlib.readaligner import ReadAligner
from reademptionlib.readalignerstats import ReadAlignerStats
from reademptionlib.readalignerstatstable import ReadAlignerStatsTable
from reademptionlib.readprocessor import ReadProcessor
from reademptionlib.sambamconverter import SamToBamConverter
from reademptionlib.wiggle import WiggleWriter


class Controller(object):

    """Manage the actions of the subcommands.

    The Controller takes care of providing the arguments like path
    names and the parallel processing of tasks.

    """

    def __init__(self, args):
        """Create an instance."""
        self._args = args
        self._paths = Paths(args.project_path)
        self._read_files = None
        self._ref_seq_files = None

    def create_project(self, version):
        """Create a new project."""
        sys.stdout.write(
            "   ___  _______   ___                 __  _\n"
            "  / _ \\/ __/ _ | / _ \\___ __ _  ___  / /_(_)__  ___\n"
            " / , _/ _// __ |/ // / -_)  ' \\/ _ \\/ __/ / _ \\/ _ \\\n"
            "/_/|_/___/_/ |_/____/\\__/_/_/_/ .__/\\__/_/\\___/_//_/\n"
            "                             / /\n"
            "====================================================\n"
            "========================================\n"
            "=======================\n"
            "==============\n\n"
            "[https://reademption.readthedocs.io/en/latest/]\n\n"
        )
        project_creator = ProjectCreator()
        project_creator.create_root_folder(self._args.project_path)
        project_creator.create_subfolders(self._paths.required_folders())
        project_creator.create_version_file(self._paths.version_path, version)
        sys.stdout.write(
            'Created folder "%s" and required subfolders.\n'
            % (self._args.project_path)
        )
        sys.stdout.write(
            'Please copy read files into folder "%s" and '
            'reference sequences files into folder "%s".\n'
            % (self._paths.read_fasta_folder, self._paths.ref_seq_folder)
        )

    def align_reads(self):
        """Perform the alignment of the reads."""
        self._args.realign = False
        self._test_folder_existance(
            self._paths.required_read_alignment_folders()
        )
        assert self._args.paired_end in [True, False]
        self._ref_seq_files = self._paths.get_ref_seq_files()
        self._paths.set_ref_seq_paths(self._ref_seq_files)
        self._test_align_file_existance()
        if not self._args.paired_end:
            # Single end reads
            self._read_files = self._paths.get_read_files()
            self._lib_names = self._paths.get_lib_names_single_end()
            self._paths.set_read_files_dep_file_lists_single_end(
                self._read_files, self._lib_names
            )
            if not self._args.realign:
                self._set_primary_aligner_paths_to_final_paths()
            self._prepare_reads_single_end()
            self._align_single_end_reads()
        else:
            # Paired end reads
            self._read_file_pairs = self._paths.get_read_file_pairs()
            self._lib_names = self._paths.get_lib_names_paired_end()
            self._paths.set_read_files_dep_file_lists_paired_end(
                self._read_file_pairs, self._lib_names
            )
            if not self._args.realign:
                self._set_primary_aligner_paths_to_final_paths()
            self._prepare_reads_paired_end()
            self._align_paired_end_reads()
        # self._sam_to_bam(
        #    self._paths.primary_read_aligner_sam_paths,
        #    self._paths.primary_read_aligner_bam_prefix_paths,
        #    self._paths.primary_read_aligner_bam_paths)
        self._generate_read_alignment_stats(
            self._lib_names,
            self._paths.primary_read_aligner_bam_paths,
            self._paths.unaligned_reads_paths,
            self._paths.primary_read_aligner_stats_path,
        )
        final_unaligned_reads_paths = self._paths.unaligned_reads_paths
        if self._args.realign:
            self._run_realigner_and_process_alignments()
            self._merge_bam_files()
            final_unaligned_reads_paths = (
                self._paths.realigned_unaligned_reads_paths
            )
        if self._args.crossalign_cleaning_str is not None:
            self._remove_crossaligned_reads()
        self._generate_read_alignment_stats(
            self._lib_names,
            self._paths.read_alignment_bam_paths,
            final_unaligned_reads_paths,
            self._paths.read_alignments_stats_path,
        )
        self._write_alignment_stat_table()

    def _remove_crossaligned_reads(self):
        self._string_to_species_and_sequence_ids()
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes
        ) as executor:
            for (
                bam_path,
                bam_with_crossmappings_path,
                bam_cleaned_tmp_path,
                crossmapped_reads_path,
            ) in zip(
                self._paths.read_alignment_bam_paths,
                self._paths.read_alignment_bam_with_crossmappings_paths,
                self._paths.read_alignment_bam_cross_cleaned_tmp_paths,
                self._paths.crossmapped_reads_paths,
            ):
                jobs.append(
                    executor.submit(
                        self._remove_crossaligned_reads_for_lib,
                        bam_path,
                        bam_with_crossmappings_path,
                        bam_cleaned_tmp_path,
                        crossmapped_reads_path,
                    )
                )
        # Evaluate thread outcome
        self._check_job_completeness(jobs)

    def _remove_crossaligned_reads_for_lib(
        self,
        bam_path,
        bam_with_crossmappings_path,
        bam_cleaned_tmp_path,
        crossmapped_reads_path,
    ):
        # Perform the removal or cross aligned reads
        cross_align_filter = CrossAlignFilter(
            bam_path,
            bam_cleaned_tmp_path,
            crossmapped_reads_path,
            self._species_and_sequence_ids,
        )
        cross_align_filter.determine_crossmapped_reads()
        cross_align_filter.write_crossmapping_free_bam()
        # Rename the original mapping file that potentially
        # contains cross aligned reads
        os.rename(bam_path, bam_with_crossmappings_path)
        os.rename(bam_path + ".bai", bam_with_crossmappings_path + ".bai")
        # Move the cross aligned filtered file to the final mapping
        # path
        os.rename(bam_cleaned_tmp_path, bam_path)
        os.rename(bam_cleaned_tmp_path + ".bai", bam_path + ".bai")

    def _string_to_species_and_sequence_ids(self):
        self._species_and_sequence_ids = {}
        orgs_and_seq_ids_strs = self._args.crossalign_cleaning_str.split(";")
        if len(orgs_and_seq_ids_strs) < 2:
            self._write_err_msg_and_quit(
                "Error! Only one organism is defined for the cross align "
                "removal. This does not make sense.\nYou gave the "
                "following input:\n%s\n" % self._args.crossalign_cleaning_str
            )
        for org_and_seq_ids_str in orgs_and_seq_ids_strs:
            org, seq_ids_str = org_and_seq_ids_str.strip().split(":")
            seq_ids = [seq_id.strip() for seq_id in seq_ids_str.split(",")]
            if "" in seq_ids:
                seq_ids.remove("")
            if len(seq_ids) < 1:
                self._write_err_msg_and_quit(
                    "Error! No sequence ID was given for the species '%s'. "
                    "This does not make sense.\nYou gave the "
                    "following input:\n%s\n"
                    % (org, self._args.crossalign_cleaning_str)
                )
            self._species_and_sequence_ids[org] = seq_ids

    def _set_primary_aligner_paths_to_final_paths(self):
        # If no remapping is performed the paths of the final bam files
        # is the paths of the primary mapper
        self._paths.primary_read_aligner_bam_prefix_paths = (
            self._paths.read_alignment_bam_prefix_paths
        )
        self._paths.primary_read_aligner_bam_paths = (
            self._paths.read_alignment_bam_paths
        )
        self._paths.primary_read_aligner_stats_path = (
            self._paths.read_alignments_stats_path
        )

    def _merge_bam_files(self):
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes
        ) as executor:
            for merged_bam, primary_aligner_bam, realigner_bam in zip(
                self._paths.read_alignment_bam_paths,
                self._paths.primary_read_aligner_bam_paths,
                self._paths.read_realigner_bam_paths,
            ):
                bam_merger = BamMerger()
                jobs.append(
                    executor.submit(
                        bam_merger.merge,
                        merged_bam,
                        primary_aligner_bam,
                        realigner_bam,
                    )
                )
        self._check_job_completeness(jobs)
        if not self._args.keep_original_alignments:
            for bam_file_list in [
                self._paths.primary_read_aligner_bam_paths,
                self._paths.read_realigner_bam_paths,
            ]:
                for bam_file in bam_file_list:
                    os.remove(bam_file)
                    os.remove("%s.bai" % bam_file)

    def _run_realigner_and_process_alignments(self):
        # As the realigner needs a *sorted* SAM file
        self._generate_sorted_tmp_sam_file()
        self._realign_unmapped_reads()
        # self._sam_to_bam(
        #    self._paths.read_realigner_sam_paths,
        #    self._paths.read_realigner_bam_prefixes_paths,
        #    self._paths.read_realigner_sam_paths)
        self._generate_read_alignment_stats(
            self._lib_names,
            self._paths.read_realigner_bam_paths,
            self._paths.realigned_unaligned_reads_paths,
            self._paths.read_realigner_stats_path,
        )

    def _test_align_file_existance(self):
        """Test if the input file for the the align subcommand exist."""
        if len(self._paths.get_read_files()) == 0:
            self._write_err_msg_and_quit("Error! No read libraries given!\n")
        if len(self._ref_seq_files) == 0:
            self._write_err_msg_and_quit(
                "Error! No reference sequence files given!\n"
            )

    def _test_folder_existance(self, task_specific_folders):
        """Test the existance of required folders."""
        for folder in (
            self._paths.required_base_folders() + task_specific_folders
        ):
            if not os.path.exists(folder):
                self._write_err_msg_and_quit(
                    "Error! Folder '%s' does not exist! Is the given project "
                    "folder name correct?\n" % folder
                )

    def _file_needs_to_be_created(self, file_path, quiet=False):
        """Test if a file exists of need to be created."""
        if not self._args.check_for_existing_files:
            return True
        if os.path.exists(file_path):
            if not quiet:
                sys.stderr.write(
                    "File %s exists. Skipping its generation.\n" % file_path
                )
            return False
        return True

    def _prepare_reads_single_end(self):
        """Manage the prepartion of reads before the actual mappings."""
        read_files_and_jobs = {}
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes
        ) as executor:
            for lib_name, read_path, processed_read_path in zip(
                self._lib_names,
                self._paths.read_paths,
                self._paths.processed_read_paths,
            ):
                if not self._file_needs_to_be_created(processed_read_path):
                    continue
                read_processor = ReadProcessor(
                    poly_a_clipping=self._args.poly_a_clipping,
                    min_read_length=self._args.min_read_length,
                    fastq=self._args.fastq,
                    min_phred_score=self._args.min_phred_score,
                    adapter=self._args.adapter,
                    reverse_complement=self._args.reverse_complement,
                )
                read_files_and_jobs[lib_name] = executor.submit(
                    read_processor.process_single_end,
                    read_path,
                    processed_read_path,
                )
        self._evaluet_job_and_generate_stat_file(read_files_and_jobs)

    def _prepare_reads_paired_end(self):
        read_files_and_jobs = {}
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes
        ) as executor:
            for lib_name, read_path_pair, processed_read_path_pair in zip(
                self._lib_names,
                self._paths.read_path_pairs,
                self._paths.processed_read_path_pairs,
            ):
                for processed_read_path in processed_read_path_pair:
                    if not self._file_needs_to_be_created(processed_read_path):
                        continue
                    read_processor = ReadProcessor(
                        poly_a_clipping=False,
                        min_read_length=self._args.min_read_length,
                        fastq=self._args.fastq,
                        min_phred_score=self._args.min_phred_score,
                        adapter=self._args.adapter,
                        reverse_complement=self._args.reverse_complement,
                    )
                    read_files_and_jobs[lib_name] = executor.submit(
                        read_processor.process_paired_end,
                        read_path_pair,
                        processed_read_path_pair,
                    )
        self._evaluet_job_and_generate_stat_file(read_files_and_jobs)

    def _evaluet_job_and_generate_stat_file(self, read_files_and_jobs):
        raw_stat_data_writer = RawStatDataWriter(pretty=True)
        # Evaluate thread outcome
        self._check_job_completeness(read_files_and_jobs.values())
        if not self._file_needs_to_be_created(
            self._paths.read_processing_stats_path
        ):
            return
        # Create a dict of the read file names and the processing
        # counting results
        read_files_and_stats = dict(
            [
                (lib_name, job.result())
                for lib_name, job in read_files_and_jobs.items()
            ]
        )
        raw_stat_data_writer.write(
            read_files_and_stats, self._paths.read_processing_stats_path
        )

    def _align_single_end_reads(self):
        """Manage the actual alignment of single end reads."""
        read_aligner = ReadAligner(self._args.segemehl_bin, self._args.progress)
        if self._file_needs_to_be_created(self._paths.index_path):
            read_aligner.build_index(
                self._paths.ref_seq_paths, self._paths.index_path
            )
        for read_path, output_path, nomatch_path, bam_path in zip(
            self._paths.processed_read_paths,
            self._paths.primary_read_aligner_bam_paths,
            self._paths.unaligned_reads_paths,
            self._paths.read_alignment_bam_paths,
        ):
            if not self._file_needs_to_be_created(output_path):
                continue
            elif not self._file_needs_to_be_created(bam_path):
                continue
            read_aligner.run_alignment(
                read_path,
                self._paths.index_path,
                self._paths.ref_seq_paths,
                output_path,
                nomatch_path,
                int(self._args.processes),
                int(self._args.segemehl_accuracy),
                float(self._args.segemehl_evalue),
                self._args.split,
                paired_end=False,
            )

    def _align_paired_end_reads(self):
        """Manage the actual alignemnt of paired end reads."""
        read_aligner = ReadAligner(self._args.segemehl_bin, self._args.progress)
        if self._file_needs_to_be_created(self._paths.index_path):
            read_aligner.build_index(
                self._paths.ref_seq_paths, self._paths.index_path
            )
        for read_path_pair, output_path, nomatch_path in zip(
            self._paths.processed_read_path_pairs,
            self._paths.primary_read_aligner_bam_paths,
            self._paths.unaligned_reads_paths,
        ):
            if not self._file_needs_to_be_created(output_path):
                continue
            # elif not self._file_needs_to_be_created(bam_path):
            #    continue
            read_aligner.run_alignment(
                read_path_pair,
                self._paths.index_path,
                self._paths.ref_seq_paths,
                output_path,
                nomatch_path,
                int(self._args.processes),
                int(self._args.segemehl_accuracy),
                float(self._args.segemehl_evalue),
                self._args.split,
                paired_end=True,
            )

    def _realign_unmapped_reads(self):
        read_realigner = ReadRealigner(self._args.lack_bin, self._args.progress)
        for (
            query_fasta_path,
            query_sam_path,
            realignment_sam_path,
            unaligned_reads_path,
        ) in zip(
            self._paths.unaligned_reads_paths,
            self._paths.read_realigner_tmp_sam_paths,
            self._paths.read_realigner_sam_paths,
            self._paths.realigned_unaligned_reads_paths,
        ):
            read_realigner.run_alignment(
                query_fasta_path,
                query_sam_path,
                self._paths.ref_seq_paths,
                realignment_sam_path,
                unaligned_reads_path,
                int(self._args.processes),
                int(self._args.segemehl_accuracy),
            )
            os.remove(query_sam_path)

    def _sam_to_bam(self, sam_paths, bam_prefixes_paths, bam_paths):
        """Manage the conversion of mapped read from SAM to BAM format."""
        sam_to_bam_converter = SamToBamConverter()
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes
        ) as executor:
            for sam_path, bam_prefix_path, bam_path in zip(
                sam_paths, bam_prefixes_paths, bam_paths
            ):
                if not self._file_needs_to_be_created(bam_path):
                    continue
                jobs.append(
                    executor.submit(
                        sam_to_bam_converter.sam_to_bam,
                        sam_path,
                        bam_prefix_path,
                    )
                )
        # Evaluate thread outcome
        self._check_job_completeness(jobs)

    def _generate_sorted_tmp_sam_file(self):
        sam_to_bam_converter = SamToBamConverter()
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes
        ) as executor:
            for bam_path, sam_path in zip(
                self._paths.primary_read_aligner_bam_paths,
                self._paths.read_realigner_tmp_sam_paths,
            ):
                jobs.append(
                    executor.submit(
                        sam_to_bam_converter.bam_to_sam, bam_path, sam_path
                    )
                )
        # Evaluate thread outcome
        self._check_job_completeness(jobs)

    def _generate_read_alignment_stats(
        self,
        lib_names,
        result_bam_paths,
        unaligned_reads_paths,
        output_stats_path,
    ):
        """Manage the generation of alingment statistics."""
        raw_stat_data_writer = RawStatDataWriter(pretty=True)
        read_files_and_jobs = {}
        if not self._file_needs_to_be_created(output_stats_path):
            return
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes
        ) as executor:
            for (
                lib_name,
                read_alignment_bam_path,
                unaligned_reads_path,
            ) in zip(lib_names, result_bam_paths, unaligned_reads_paths):
                read_aligner_stats = ReadAlignerStats()
                read_files_and_jobs[lib_name] = executor.submit(
                    read_aligner_stats.count,
                    read_alignment_bam_path,
                    unaligned_reads_path,
                )
        # Evaluate thread outcome
        self._check_job_completeness(read_files_and_jobs.values())
        read_files_and_stats = dict(
            [
                (lib_name, job.result())
                for lib_name, job in read_files_and_jobs.items()
            ]
        )
        raw_stat_data_writer.write(read_files_and_stats, output_stats_path)

    def _write_alignment_stat_table(self):
        """Manage the creation of the mapping statistic output table."""
        raw_stat_data_reader = RawStatDataReader()
        read_processing_stats = raw_stat_data_reader.read(
            self._paths.read_processing_stats_path
        )
        final_alignment_stats = raw_stat_data_reader.read(
            self._paths.read_alignments_stats_path
        )
        realignment_stats = None
        primary_aligner_stats = None
        if self._args.realign:
            primary_aligner_stats = raw_stat_data_reader.read(
                self._paths.primary_read_aligner_stats_path
            )
            realignment_stats = raw_stat_data_reader.read(
                self._paths.read_realigner_stats_path
            )
        read_aligner_stats_table = ReadAlignerStatsTable(
            read_processing_stats,
            final_alignment_stats,
            primary_aligner_stats,
            realignment_stats,
            self._lib_names,
            self._paths.read_alignment_stats_table_path,
            self._args.paired_end,
        )
        read_aligner_stats_table.write()

    def _ref_ids_to_file(self, ref_seq_paths):
        """Translate the reference ID to file paths."""
        ref_ids_to_file = {}
        fasta_parser = FastaParser()
        for ref_seq_path in ref_seq_paths:
            ref_seq_file = os.path.basename(ref_seq_path)
            with open(ref_seq_path) as ref_seq_fh:
                ref_seq_id = fasta_parser.header_id(
                    fasta_parser.single_entry_file_header(ref_seq_fh)
                )
                ref_ids_to_file[ref_seq_id] = ref_seq_file
        return ref_ids_to_file

    def create_coverage_files(self):
        """Create coverage files based on the read alignments.

        The coverages are calculated per replicon and the results are
        written to the output file. This might be slower but if all
        coverages are determined at once the data structure will become
        too large when working with large reference sequences.

        """
        self._test_folder_existance(self._paths.required_coverage_folders())
        raw_stat_data_reader = RawStatDataReader()
        alignment_stats = [
            raw_stat_data_reader.read(self._paths.read_alignments_stats_path)
        ]
        lib_names = list(alignment_stats[0].keys())
        was_paired_end_alignment = self._was_paired_end_alignment(lib_names)
        if not was_paired_end_alignment:
            self._paths.set_read_files_dep_file_lists_single_end(
                self._paths.get_read_files(), lib_names
            )
        else:
            self._paths.set_read_files_dep_file_lists_paired_end(
                self._paths.get_read_files(), lib_names
            )
        # Get number of aligned or number of uniquely aligned reads
        if not self._args.normalize_by_uniquely:
            aligned_counting = "no_of_aligned_reads"
        else:
            aligned_counting = "no_of_uniquely_aligned_reads"
        read_files_aligned_read_freq = dict(
            [
                (read_file, round(attributes["stats_total"][aligned_counting]))
                for read_file, attributes in alignment_stats[0].items()
            ]
        )
        min_no_of_aligned_reads = float(
            min(read_files_aligned_read_freq.values())
        )
        # Run the generation of coverage in parallel
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes
        ) as executor:
            for lib_name, bam_path in zip(
                lib_names, self._paths.read_alignment_bam_paths
            ):
                no_of_aligned_reads = float(
                    read_files_aligned_read_freq[lib_name]
                )
                jobs.append(
                    executor.submit(
                        self._create_coverage_files_for_lib,
                        lib_name,
                        bam_path,
                        no_of_aligned_reads,
                        min_no_of_aligned_reads,
                    )
                )
        # Evaluate thread outcome
        self._check_job_completeness(jobs)

    def _all_coverage_file_exist(
        self, lib_name, strands, no_of_aligned_reads, min_no_of_aligned_reads
    ):
        """Test the existance of all coverage files of a library"""
        files = []
        for strand in strands:
            files.append(self._paths.wiggle_file_raw_path(lib_name, strand))
            files.append(
                self._paths.wiggle_file_tnoar_norm_min_path(
                    lib_name,
                    strand,
                    multi=min_no_of_aligned_reads,
                    div=no_of_aligned_reads,
                )
            )
            files.append(
                self._paths.wiggle_file_tnoar_norm_mil_path(
                    lib_name, strand, multi=1000000, div=no_of_aligned_reads
                )
            )
        if not any(
            [self._file_needs_to_be_created(file, quiet=True) for file in files]
        ):
            sys.stderr.write(
                "The files %s exists. Skipping their generation.\n"
                % ", ".join(files)
            )
            return True
        return False

    def _create_coverage_files_for_lib(
        self, lib_name, bam_path, no_of_aligned_reads, min_no_of_aligned_reads
    ):
        """Perform the coverage calculation for a given library."""
        if not self._args.non_strand_specific:
            strands = ["forward", "reverse"]
        else:
            strands = ["forward_and_reverse"]
        if self._all_coverage_file_exist(
            lib_name, strands, no_of_aligned_reads, min_no_of_aligned_reads
        ):
            return
        read_count_splitting = True
        if self._args.skip_read_count_splitting:
            read_count_splitting = False
        coverage_calculator = CoverageCalculator(
            read_count_splitting=read_count_splitting,
            uniquely_aligned_only=self._args.unique_only,
            coverage_style=self._args.coverage_style,
            clip_length=self._args.clip_length,
            non_strand_specific=self._args.non_strand_specific,
        )
        (
            coverage_writers_raw,
            coverage_writers_tnoar_min_norm,
            coverage_writers_tnoar_mil_norm,
        ) = self._wiggle_writers(
            lib_name, strands, no_of_aligned_reads, min_no_of_aligned_reads
        )
        for ref_seq, coverages in coverage_calculator.ref_seq_and_coverages(
            bam_path
        ):
            for strand in strands:
                coverage_writers_raw[strand].write_replicons_coverages(
                    ref_seq, coverages[strand]
                )
                coverage_writers_tnoar_min_norm[
                    strand
                ].write_replicons_coverages(
                    ref_seq,
                    coverages[strand],
                    factor=min_no_of_aligned_reads / no_of_aligned_reads,
                )
                coverage_writers_tnoar_mil_norm[
                    strand
                ].write_replicons_coverages(
                    ref_seq,
                    coverages[strand],
                    factor=1000000 / no_of_aligned_reads,
                )
        for strand in strands:
            coverage_writers_raw[strand].close_file()

    def _wiggle_writers(
        self, lib_name, strands, no_of_aligned_reads, min_no_of_aligned_reads
    ):
        """Write the calculated coverages to wiggle files."""
        coverage_writers_raw = dict(
            [
                (
                    strand,
                    WiggleWriter(
                        "%s_%s" % (lib_name, strand),
                        open(
                            self._paths.wiggle_file_raw_path(lib_name, strand),
                            "w",
                        ),
                    ),
                )
                for strand in strands
            ]
        )
        coverage_writers_tnoar_min_norm = dict(
            [
                (
                    strand,
                    WiggleWriter(
                        "%s_%s" % (lib_name, strand),
                        open(
                            self._paths.wiggle_file_tnoar_norm_min_path(
                                lib_name,
                                strand,
                                multi=min_no_of_aligned_reads,
                                div=no_of_aligned_reads,
                            ),
                            "w",
                        ),
                    ),
                )
                for strand in strands
            ]
        )
        coverage_writers_tnoar_mil_norm = dict(
            [
                (
                    strand,
                    WiggleWriter(
                        "%s_%s" % (lib_name, strand),
                        open(
                            self._paths.wiggle_file_tnoar_norm_mil_path(
                                lib_name,
                                strand,
                                multi=1000000,
                                div=no_of_aligned_reads,
                            ),
                            "w",
                        ),
                    ),
                )
                for strand in strands
            ]
        )
        return (
            coverage_writers_raw,
            coverage_writers_tnoar_min_norm,
            coverage_writers_tnoar_mil_norm,
        )

    def _check_job_completeness(self, jobs):
        """Check the completness of each job in a list"""
        for job in concurrent.futures.as_completed(jobs):
            if job.exception():
                raise (job.exception())

    def quantify_gene_wise(self):
        """Manage the counting of aligned reads per gene."""
        self._test_folder_existance(self._paths.required_gene_quanti_folders())
        norm_by_alignment_freq = True
        norm_by_overlap_freq = True
        if self._args.no_count_split_by_alignment_no:
            norm_by_alignment_freq = False
        if self._args.no_count_splitting_by_gene_no:
            norm_by_overlap_freq = False
        raw_stat_data_reader = RawStatDataReader()
        alignment_stats = [
            raw_stat_data_reader.read(self._paths.read_alignments_stats_path)
        ]
        lib_names = sorted(list(alignment_stats[0].keys()))
        annotation_files = self._paths.get_annotation_files()
        self._paths.set_annotation_paths(annotation_files)
        was_paired_end_alignment = self._was_paired_end_alignment(lib_names)
        if not was_paired_end_alignment:
            self._paths.set_read_files_dep_file_lists_single_end(
                self._paths.get_read_files(), lib_names
            )
        else:
            self._paths.set_read_files_dep_file_lists_paired_end(
                self._paths.get_read_files(), lib_names
            )
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes
        ) as executor:
            for lib_name, read_alignment_path in zip(
                lib_names, self._paths.read_alignment_bam_paths
            ):
                jobs.append(
                    executor.submit(
                        self._quantify_gene_wise,
                        lib_name,
                        read_alignment_path,
                        norm_by_alignment_freq,
                        norm_by_overlap_freq,
                        annotation_files,
                    )
                )
        # Evaluate thread outcome
        self._check_job_completeness(jobs)
        self._gene_quanti_create_overview(
            annotation_files, self._paths.annotation_paths, lib_names
        )

    def _was_paired_end_alignment(self, lib_names):
        """Check if the mapping was done in paired- or single-end mode"""
        if len(lib_names) * 2 == len(self._paths.get_read_files()):
            return True
        return False

    def _quantify_gene_wise(
        self,
        lib_name,
        read_alignment_path,
        norm_by_alignment_freq,
        norm_by_overlap_freq,
        annotation_files,
    ):
        """Perform the gene wise quantification for a given library."""
        gene_quanti_paths = [
            self._paths.gene_quanti_path(lib_name, annotation_file)
            for annotation_file in annotation_files
        ]
        # Check if all output files for this library exist - if so
        # skip their creation
        if not any(
            [
                self._file_needs_to_be_created(gene_quanti_path, quiet=True)
                for gene_quanti_path in gene_quanti_paths
            ]
        ):
            sys.stderr.write(
                "The file(s) %s exist(s). Skipping their/its generation.\n"
                % ", ".join(gene_quanti_paths)
            )
            return
        strand_specific = True
        if self._args.non_strand_specific:
            strand_specific = False
        gene_wise_quantification = GeneWiseQuantification(
            min_overlap=self._args.min_overlap,
            read_region=self._args.read_region,
            clip_length=self._args.clip_length,
            norm_by_alignment_freq=norm_by_alignment_freq,
            norm_by_overlap_freq=norm_by_overlap_freq,
            allowed_features_str=self._args.allowed_features,
            add_antisense=self._args.add_antisense,
            antisense_only=self._args.antisense_only,
            strand_specific=strand_specific,
            unique_only=self._args.unique_only,
        )
        if norm_by_overlap_freq:
            gene_wise_quantification.calc_overlaps_per_alignment(
                read_alignment_path, self._paths.annotation_paths
            )
        for annotation_file, annotation_path in zip(
            annotation_files, self._paths.annotation_paths
        ):
            gene_wise_quantification.quantify(
                read_alignment_path,
                annotation_path,
                self._paths.gene_quanti_path(lib_name, annotation_file),
                self._args.pseudocounts,
            )

    def _gene_quanti_create_overview(
        self, annotation_files, annotation_paths, lib_names
    ):
        """Create an overview table of all gene quantification for all libs."""
        strand_specific = True
        if self._args.non_strand_specific:
            strand_specific = False
        gene_wise_overview = GeneWiseOverview(
            allowed_features_str=self._args.allowed_features,
            add_antisense=self._args.add_antisense,
            antisense_only=self._args.antisense_only,
            strand_specific=strand_specific,
        )
        path_and_name_combos = {}
        for annotation_file, annotation_path in zip(
            annotation_files, annotation_paths
        ):
            path_and_name_combos[annotation_path] = []
            for read_file in lib_names:
                path_and_name_combos[annotation_path].append(
                    [
                        read_file,
                        self._paths.gene_quanti_path(
                            read_file, annotation_file
                        ),
                    ]
                )
        if self._file_needs_to_be_created(
            self._paths.gene_wise_quanti_combined_path
        ):
            gene_wise_overview.create_overview_raw_countings(
                path_and_name_combos,
                lib_names,
                self._paths.gene_wise_quanti_combined_path,
            )
        if self._file_needs_to_be_created(
            self._paths.gene_wise_quanti_combined_rpkm_path
        ):
            gene_wise_overview.create_overview_rpkm(
                path_and_name_combos,
                lib_names,
                self._paths.gene_wise_quanti_combined_rpkm_path,
                self._libs_and_total_num_of_aligned_reads(),
            )
        if self._file_needs_to_be_created(
            self._paths.gene_wise_quanti_combined_tnoar_path
        ):
            gene_wise_overview.create_overview_norm_by_tnoar(
                path_and_name_combos,
                lib_names,
                self._paths.gene_wise_quanti_combined_tnoar_path,
                self._libs_and_total_num_of_aligned_reads(),
            )
        if self._file_needs_to_be_created(
            self._paths.gene_wise_quanti_combined_tpm_path
        ):
            gene_wise_overview.create_overview_tpm(
                self._paths.gene_wise_quanti_combined_path,
                self._paths.gene_wise_quanti_combined_tpm_path,
            )

    def _libs_and_total_num_of_aligned_reads(self):
        """Read the total number of reads per library."""
        with open(
            self._paths.read_alignments_stats_path
        ) as read_aligner_stats_fh:
            read_aligner_stats = json.loads(read_aligner_stats_fh.read())
        return dict(
            [
                (lib, values["stats_total"]["no_of_aligned_reads"])
                for lib, values in read_aligner_stats.items()
            ]
        )

    def _libs_and_total_num_of_uniquely_aligned_reads(self):
        """Read the total number of reads per library."""
        with open(
            self._paths.read_alignments_stats_path
        ) as read_aligner_stats_fh:
            read_aligner_stats = json.loads(read_aligner_stats_fh.read())
        return dict(
            [
                (lib, values["stats_total"]["no_of_uniquely_aligned_reads"])
                for lib, values in read_aligner_stats.items()
            ]
        )

    def compare_with_deseq(self):
        """Manage the pairwise expression comparison with DESeq."""
        self._test_folder_existance(self._paths.required_deseq_folders())
        arg_libs = [
            self._paths._clean_file_name(lib)
            for lib in self._args.libs.split(",")
        ]
        conditions = self._args.conditions.split(",")
        self._check_deseq_args(arg_libs, conditions)
        deseq_runner = DESeqRunner(
            arg_libs,
            conditions,
            self._paths.deseq_raw_folder,
            self._paths.deseq_extended_folder,
            self._paths.deseq_script_path,
            self._paths.deseq_pca_heatmap_path,
            self._paths.gene_wise_quanti_combined_path,
            self._paths.deseq_tmp_session_info_script,
            self._paths.deseq_session_info,
            self._args.fc_shrinkage_off,
            self._args.cooks_cutoff_off,
        )
        deseq_runner.create_deseq_script_file()
        deseq_runner.write_session_info_file()
        deseq_runner.run_deseq()
        deseq_runner.merge_counting_files_with_results()

    def _check_deseq_args(self, arg_libs, conditions):
        """Test if the given arguments are sufficient."""
        if len(arg_libs) != len(conditions):
            self._write_err_msg_and_quit(
                "Error - The read library file list and condition list must "
                "have the same number of elements. You entered \n%s "
                "(= %s elements)\nand \n%s (= %s elements).\n"
                % (
                    self._args.libs,
                    len(arg_libs),
                    self._args.conditions,
                    len(conditions),
                )
            )
        raw_stat_data_reader = RawStatDataReader()
        alignment_stats = [
            raw_stat_data_reader.read(self._paths.read_alignments_stats_path)
        ]
        lib_names = list(alignment_stats[0].keys())
        if len(lib_names) != len(arg_libs):
            self._write_err_msg_and_quit(
                "The number of read libraries is lower or higher than "
                "expected. The following read libs are available: %s\nThe "
                'following read list string is suggested: "%s"\n'
                % (", ".join(lib_names), ",".join(lib_names))
            )
        for lib in lib_names:
            if lib not in arg_libs:
                self._write_err_msg_and_quit(
                    'The library "%s" is not present in your list of '
                    "libraries. Please add it.\n" % (lib)
                )

    def _write_err_msg_and_quit(self, msg):
        """Write error message and close the program gracefully."""
        sys.stderr.write(msg)
        sys.exit(1)

    def viz_align(self):
        """Generate plots based on the read processing and mapping"""
        from reademptionlib.vizalign import AlignViz

        align_viz = AlignViz(
            self._paths.get_lib_names_single_end()
            if not self._args.paired_end
            else self._paths.get_lib_names_paired_end(),
            self._paths.read_processing_stats_path,
            self._paths.read_alignments_stats_path,
        )
        align_viz.read_stat_files()
        align_viz.plot_input_read_length(
            self._paths.viz_align_input_read_length_plot_path
        )
        align_viz.plot_processed_read_length(
            self._paths.viz_align_processed_reads_length_plot_path
        )

    def viz_gene_quanti(self):
        """Generate plots based on the gene-wise read countings"""
        from reademptionlib.vizgenequanti import GeneQuantiViz

        gene_quanti_viz = GeneQuantiViz(
            self._paths.gene_wise_quanti_combined_path,
            self._paths.get_lib_names_single_end()
            if not self._args.paired_end
            else self._paths.get_lib_names_paired_end(),
        )
        gene_quanti_viz.parse_input_table()
        gene_quanti_viz.plot_correlations(
            self._paths.viz_gene_quanti_scatter_plot_path
        )
        gene_quanti_viz.plot_annotation_class_quantification(
            self._paths.viz_gene_quanti_rna_classes_plot_path
        )

    def viz_deseq(self):
        """Generate plots based on the DESeq analysis"""
        from reademptionlib.vizdeseq import DESeqViz

        deseq_path_template = (
            self._paths.deseq_raw_folder + "/deseq_comp_%s_vs_%s.csv"
        )
        deseq_viz = DESeqViz(
            self._paths.deseq_script_path,
            deseq_path_template,
            max_pvalue=self._args.max_pvalue,
        )
        deseq_viz.create_scatter_plots(self._paths.viz_deseq_scatter_plot_path)
        deseq_viz.create_volcano_plots(
            self._paths.viz_deseq_volcano_plot_path,
            self._paths.viz_deseq_volcano_plot_adj_path,
        )
