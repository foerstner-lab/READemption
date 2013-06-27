import concurrent.futures
import os
import sys
import json
sys.path.append(".")
from libs.coveragecalculator import CoverageCalculator
from libs.deseq import DESeqRunner
from libs.fasta import FastaParser
from libs.genewisequanti import GeneWiseQuantification, GeneWiseOverview
from libs.parameterlog import ParameterLogger
from libs.paths import Paths
from libs.projectcreator import ProjectCreator
from libs.rawstatdata import RawStatDataWriter, RawStatDataReader
from libs.readaligner import ReadAligner
from libs.readalignerstats import ReadAlignerStats
from libs.readalignerstatstable import ReadAlignerStatsTable
from libs.readprocessor import ReadProcessor
from libs.sambamconverter import SamToBamConverter
from libs.wiggle import WiggleWriter

class Controller(object):

    def __init__(self, args):
        """Create an instance."""
        self.args = args
        self.paths = Paths(args.project_path)

    def create_project(self, version):
        """Create a new project."""
        project_creator = ProjectCreator()
        project_creator.create_root_folder(self.args.project_path)
        project_creator.create_subfolders(self.paths.required_folders())
        project_creator.create_version_file(self.paths.version_path, version)
        sys.stdout.write("Created folder \"%s\" and required subfolders.\n" % (
                self.args.project_path))
        sys.stdout.write("Please copy read files into folder \"%s\" and "
                         "reference sequences files into folder \"%s\".\n" % (
                self.paths.read_fasta_folder, self.paths.ref_seq_folder))

    def align_reads(self):
        """Perform the alignment of the reads."""
        self._test_folder_existance(
            self.paths.required_read_alignment_folders())
        self.read_files = self.paths.get_read_files()
        self.lib_names = self.paths.get_lib_names()
        ref_seq_files = self.paths.get_ref_seq_files()
        self.paths.set_read_files_dep_file_lists(
            self.read_files, self.lib_names)
        self.paths.set_ref_seq_paths(ref_seq_files)
        self._prepare_reads()
        self._align_reads()
        self._sam_to_bam()
        self._generate_read_alignment_stats()
        self._write_alignment_stat_table()

    def _test_folder_existance(self, task_specific_folders):
        for folder in (
            self.paths.required_base_folders() + task_specific_folders):
            if not os.path.exists(folder):
                sys.stderr.write(
                    "Error! Folder '%s' does not exist! Is the given project "
                    "folder name correct?\n" % folder)
                sys.exit(2)

    def _file_needs_to_be_created(self, file_path):
        if self.args.force is True:
            return True
        if os.path.exists(file_path):
            sys.stderr.write(
                "File %s exists. Skipping its generation.\n" % file_path)
            return False
        return True

    def _prepare_reads(self):
        raw_stat_data_writer = RawStatDataWriter(pretty=True)
        read_files_and_jobs = {}
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self.args.processes) as executor:
            for lib_name, read_path, processed_read_path in zip(
                    self.lib_names, self.paths.read_paths, 
                    self.paths.processed_read_paths):
                read_processor = ReadProcessor(
                    poly_a_clipping=self.args.poly_a_clipping,
                    min_read_length=self.args.min_read_length)
                read_files_and_jobs[lib_name]  = executor.submit(
                    read_processor.process, read_path, processed_read_path)
        # Evaluate thread outcome
        self._check_job_completeness(read_files_and_jobs.values())
        # Create a dict of the read file names and the processing
        # counting results
        read_files_and_stats = dict(
            [(lib_name, job.result()) for lib_name, job in
             read_files_and_jobs.items()])
        raw_stat_data_writer.write(
            read_files_and_stats, self.paths.read_processing_stats_path)

    def _align_reads(self):
        read_aligner = ReadAligner(segemehl_bin=self.args.segemehl_bin)
        if self._file_needs_to_be_created(self.paths.index_path) is True:
            read_aligner.build_index(
                self.paths.ref_seq_paths, self.paths.index_path)
        for read_path, output_path, nomatch_path, bam_path in zip(
            self.paths.processed_read_paths, 
            self.paths.read_alignment_result_sam_paths, 
            self.paths.unaligned_reads_paths, 
            self.paths.read_alignment_result_bam_paths):
            if self._file_needs_to_be_created(output_path) is False:
                continue
            elif self._file_needs_to_be_created(bam_path) is False:
                continue
            read_aligner.run_alignment(
                read_path, self.paths.index_path, self.paths.ref_seq_paths, 
                output_path, nomatch_path, int(self.args.processes),
                int(self.args.segemehl_accuracy), 
                float(self.args.segemehl_evalue), self.args.split)

    def _sam_to_bam(self):
        sam_to_bam_converter = SamToBamConverter()
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self.args.processes) as executor:
            for sam_path, bam_prefix_path, bam_path in zip(
                    self.paths.read_alignment_result_sam_paths,
                    self.paths.read_alignment_result_bam_prefixes_paths,
                    self.paths.read_alignment_result_bam_paths):
                if self._file_needs_to_be_created(bam_path) is False:
                    continue
                jobs.append(executor.submit(
                    sam_to_bam_converter.sam_to_bam, sam_path, bam_prefix_path))
        # Evaluate thread outcome
        self._check_job_completeness(jobs)

    def _generate_read_alignment_stats(self):
        raw_stat_data_writer = RawStatDataWriter(pretty=True)
        read_files_and_jobs = {}
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self.args.processes) as executor:
            for (lib_name, read_alignment_result_bam_path,
                 unaligned_reads_path) in zip(
                     self.lib_names, 
                     self.paths.read_alignment_result_bam_paths,
                     self.paths.unaligned_reads_paths):
                read_aligner_stats = ReadAlignerStats()
                read_files_and_jobs[lib_name]  = executor.submit(
                    read_aligner_stats.count, read_alignment_result_bam_path,
                    unaligned_reads_path)
        # Evaluate thread outcome
        self._check_job_completeness(read_files_and_jobs.values())
        read_files_and_stats = dict(
            [(lib_name, job.result()) 
             for lib_name, job in read_files_and_jobs.items()])
        raw_stat_data_writer.write(
            read_files_and_stats, self.paths.read_aligner_stats_path)

    def _write_alignment_stat_table(self):
        raw_stat_data_reader = RawStatDataReader()
        read_processing_stats = raw_stat_data_reader.read(
            self.paths.read_processing_stats_path)
        alignment_stats = raw_stat_data_reader.read(
            self.paths.read_aligner_stats_path)
        read_aligner_stats_table = ReadAlignerStatsTable(
            read_processing_stats, alignment_stats, 
            self.lib_names, self.paths.read_alignment_stats_table_path)
        read_aligner_stats_table.write()

    def _ref_ids_to_file(self, ref_seq_paths):
        ref_ids_to_file = {}
        fasta_parser = FastaParser()
        for ref_seq_path in ref_seq_paths:
            ref_seq_file = os.path.basename(ref_seq_path)
            ref_seq_id = fasta_parser.header_id(
                fasta_parser.single_entry_file_header(open(ref_seq_path)))
            ref_ids_to_file[ref_seq_id] = ref_seq_file
        return ref_ids_to_file

    def create_coverage_files(self):
        """Create coverage files based on the read alignments.

        The coverages are calculated per replicon and the results are
        written to the output file. This might be slower but if all
        coveragers are detmined at once the data structure will become
        too large when working with large reference sequences.

        """
        self._test_folder_existance(
            self.paths.required_coverage_folders())
        lib_names = self.paths.get_lib_names()
        self.paths.set_read_files_dep_file_lists(
            self.paths.get_read_files(), lib_names)
        raw_stat_data_reader = RawStatDataReader()
        alignment_stats = [
            raw_stat_data_reader.read(
            self.paths.read_aligner_stats_path)]
        read_files_aligned_read_freq = dict([
            (read_file,
             round(attributes["stats_total"]["no_of_aligned_reads"]))
             for read_file, attributes in alignment_stats[0].items()])
        min_no_of_aligned_reads = float(min(
            read_files_aligned_read_freq.values()))
        # Run the generation of coverage in parallel
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self.args.processes) as executor:
            for lib_name, bam_path in zip(
                lib_names, self.paths.read_alignment_result_bam_paths):
                no_of_aligned_reads = float(
                    read_files_aligned_read_freq[lib_name])
                jobs.append(executor.submit(
                        self._create_coverage_files_for_lib,
                        lib_name, bam_path, no_of_aligned_reads,
                        min_no_of_aligned_reads))
        # Evaluate thread outcome
        self._check_job_completeness(jobs)

    def _create_coverage_files_for_lib(
            self, read_file, bam_path, no_of_aligned_reads,
            min_no_of_aligned_reads):
        strands = ["forward", "reverse"]
        read_count_splitting = True
        if self.args.skip_read_count_splitting is True:
            read_count_splitting = False
        coverage_calculator = CoverageCalculator(
            read_count_splitting=read_count_splitting,
            uniqueley_aligned_only=self.args.unique_only,
            first_base_only=self.args.first_base_only)
        (coverage_writers_raw, coverage_writers_tnoar_min_norm,
         coverage_writers_tnoar_mil_norm) = self._wiggle_writers(
             read_file, strands, no_of_aligned_reads, min_no_of_aligned_reads)
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

    def _wiggle_writers(self, read_file, strands, no_of_aligned_reads,
                        min_no_of_aligned_reads):
        coverage_writers_raw = dict([(
            strand, WiggleWriter(
                "%s_%s" % (read_file, strand),
                open(self.paths.wiggle_file_raw_path(read_file, strand), "w")))
                for strand in strands])
        coverage_writers_tnoar_min_norm = dict([(
            strand, WiggleWriter(
                "%s_%s" % (read_file, strand),
                open(self.paths.wiggle_file_tnoar_norm_min_path(
                    read_file, strand, multi=min_no_of_aligned_reads,
                    div=no_of_aligned_reads), "w")))
                for strand in strands])
        coverage_writers_tnoar_mil_norm = dict([(
            strand, WiggleWriter(
                "%s_%s" % (read_file, strand),
                open(self.paths.wiggle_file_tnoar_norm_mil_path(
                    read_file, strand, multi=1000000,
                    div=no_of_aligned_reads), "w")))
                for strand in strands])
        return (coverage_writers_raw, coverage_writers_tnoar_min_norm, 
                coverage_writers_tnoar_mil_norm)

    def _check_job_completeness(self, jobs):
        """Check the completness of each job in a list"""
        for job in concurrent.futures.as_completed(jobs):
            if job.exception():
                raise(job.exception())

    def quantify_gene_wise(self):
        self._test_folder_existance(
            self.paths.required_gene_quanti_folders())
        norm_by_alignment_freq = True
        norm_by_overlap_freq = True
        if self.args.skip_norm_by_alignment_freq:
            norm_by_alignment_freq = False
        if self.args.skip_norm_by_overlap_freq:
            norm_by_overlap_freq = False
        lib_names = self.paths.get_lib_names()
        annotation_files = self.paths.get_annotation_files()
        self.paths.set_annotation_paths(annotation_files)
        self.paths.set_read_files_dep_file_lists(
            self.paths.get_read_files(), lib_names)
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self.args.processes) as executor:
            for lib_name, read_alignment_path in zip(
                lib_names, self.paths.read_alignment_result_bam_paths):
                jobs.append(executor.submit(
                        self._quantify_gene_wise, lib_name,
                        read_alignment_path, norm_by_alignment_freq,
                        norm_by_overlap_freq, annotation_files))
        # Evaluate thread outcome
        self._check_job_completeness(jobs)
        self._gene_quanti_create_overview(
            annotation_files, self.paths.annotation_paths, lib_names)

    def _quantify_gene_wise(self, lib_name, read_alignment_path, 
                    norm_by_alignment_freq,  norm_by_overlap_freq, 
                    annotation_files):
        gene_wise_quantification = GeneWiseQuantification(
            min_overlap=self.args.min_overlap,
            norm_by_alignment_freq=norm_by_alignment_freq,
            norm_by_overlap_freq=norm_by_overlap_freq,
            allowed_features_str=self.args.allowed_features,
            skip_antisense=self.args.skip_antisense,
            unique_only=self.args.unique_only)
        gene_wise_quantification.calc_overlaps_per_alignment(
            read_alignment_path, self.paths.annotation_paths)
        for  annotation_file, annotation_path in zip(
            annotation_files, self.paths.annotation_paths):
            gene_wise_quantification.quantify(
                read_alignment_path, annotation_path,
                self.paths.gene_quanti_path(lib_name, annotation_file),
                self.args.pseudocounts)

    def _gene_quanti_create_overview(
            self, annotation_files, annotation_paths, lib_names):
        gene_wise_overview = GeneWiseOverview(
            allowed_features_str=self.args.allowed_features,
            skip_antisense=self.args.skip_antisense)
        path_and_name_combos = {}
        for annotation_file, annotation_path in zip(
                annotation_files, annotation_paths):
            path_and_name_combos[annotation_path] = []
            for read_file in lib_names:
                path_and_name_combos[annotation_path].append(
                    [read_file, self.paths.gene_quanti_path(
                        read_file, annotation_file)])
        gene_wise_overview.create_overview_raw_countings(
            path_and_name_combos, lib_names,
            self.paths.gene_wise_quanti_combined_path)
        gene_wise_overview.create_overview_rpkm(
            path_and_name_combos, lib_names,
            self.paths.gene_wise_quanti_combined_rpkm_path,
            self._libs_and_total_number_of_mapped_reads())
        gene_wise_overview.create_overview_norm_by_tnoar(
            path_and_name_combos, lib_names,
            self.paths.gene_wise_quanti_combined_tnoar_path,
            self._libs_and_total_number_of_mapped_reads())

    def _libs_and_total_number_of_mapped_reads(self):
        read_aligner_stats = json.loads(
            open(self.paths.read_aligner_stats_path).read())
        return dict([(lib, values["stats_total"]["no_of_aligned_reads"])
                     for lib, values in read_aligner_stats.items()])

    def compare_with_deseq(self):
        self._test_folder_existance(
            self.paths.required_deseq_folders())
        libs = self.args.libs.split(",")
        conditions = self.args.conditions.split(",")
        self._check_deseq_args(libs, conditions)
        deseq_runner = DESeqRunner(
            libs, conditions, self.paths.deseq_raw_folder,
            self.paths.deseq_extended_folder, self.paths.deseq_script_path,
            self.paths.gene_wise_quanti_combined_path, 
            self.paths.deseq_tmp_session_info_script,
            self.paths.deseq_session_info, no_replicates=self.args.no_replicates)
        deseq_runner.create_deseq_script_file()
        deseq_runner.write_session_info_file()
        deseq_runner.run_deseq()
        deseq_runner.merge_counting_files_with_results()

    def _check_deseq_args(self, libs, conditions):
        if len(libs) != len(conditions):
            self._write_err_msg_and_quit(
                "Error - The read library file list and condition list must "
                "have the same number of elements. You entered \n%s "
                "(= %s elements)\nand \n%s (= %s elements).\n" % (
                    self.args.libs, len(libs), self.args.conditions, len(conditions)))
        read_files = self.paths.get_read_files()
        if len(libs) != len(read_files):
            self._write_err_msg_and_quit(
                "The number of read libraries is lower or higher than "
                "expected. The following read libs are available: %s\nThe "
                "following read list string is suggested: \"%s\"\n" % (
                    ", ".join(read_files), ",".join(read_files)))
        for read_file in read_files:
            if read_file not in libs:
                self._write_err_msg_and_quit(
                    "There library \"%s\" is not given in your list of "
                    "libraries. Please add it.\n" % (read_file))

    def _write_err_msg_and_quit(self, msg):
        sys.stderr.write(msg)
        sys.exit(1)
