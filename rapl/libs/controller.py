import concurrent.futures
import os
import sys
sys.path.append(".")
from libs.fasta import FastaParser
from libs.coveragecreator import CoverageCreator
from libs.parameterlog import ParameterLogger
from libs.paths import Paths
from libs.projectcreator import ProjectCreator
from libs.readprocessor import ReadProcessor
from libs.readaligner import ReadAligner
from libs.readalignerstats import ReadAlignerStats
from libs.sambamconverter import SamToBamConverter
from libs.genewisequanti import GeneWiseQuantification, GeneWiseOverview
from libs.rawstatdata import RawStatDataWriter, RawStatDataReader

class Controller(object):

    def __init__(self, args):
        """Create an instance."""
        self.args = args
        self.paths = Paths(args.project_path)

    def create_project(self):
        """Create a new project."""
        project_creator = ProjectCreator()
        project_creator.create_root_folder(self.args.project_path)
        project_creator.create_subfolders(self.paths.required_folders())
        sys.stdout.write("Created folder \"%s\" and required subfolders.\n" % (
                self.args.project_path))
        sys.stdout.write("Please copy read files into folder \"%s\" and "
                         "reference sequences files into folder \"%s\".\n" % (
                self.paths.read_fasta_folder, self.paths.ref_seq_folder))

    def align_reads(self):
        """Perform the alignment of the reads."""
        self.read_files = self.paths._get_read_files()
        ref_seq_files = self.paths._get_ref_seq_files()
        self.paths.set_read_files_dep_file_lists(self.read_files)
        self.paths.set_ref_seq_paths(ref_seq_files)
        self._prepare_reads()
        self._align_reads()
        self._sam_to_bam()
        self._generate_read_alignment_stats()
        self._write_alignment_stat_table()

    def _prepare_reads(self):
        raw_stat_data_writer = RawStatDataWriter(pretty=True)
        read_files_and_jobs = {}
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self.args.processes) as executor:
            for read_file, read_path, processed_read_path in zip(
                    self.read_files, self.paths.read_paths,
                    self.paths.processed_read_paths):
                read_processor = ReadProcessor(
                    min_read_length=self.args.min_read_length)
                read_files_and_jobs[read_file]  = executor.submit(
                    read_processor.process, read_path,
                    processed_read_path)
        # Evaluate thread outcome
        self._check_job_completeness(read_files_and_jobs.values())
        # Create a dict of the read file names and the processing
        # counting results
        read_files_and_stats = dict(
            [(read_file, job.result()) for read_file, job in
             read_files_and_jobs.items()])
        raw_stat_data_writer.write(
            read_files_and_stats, self.paths.read_processing_stats_path)

    def _align_reads(self):
        read_aligner = ReadAligner(segemehl_bin=self.args.segemehl_bin)
        read_aligner.build_index(
            self.paths.ref_seq_paths, self.paths.index_path)
        read_aligner.run_alignment(
            self.paths.processed_read_paths,
            self.paths.ref_seq_paths, self.paths.index_path,
            self.paths.read_alignment_result_sam_paths,
            self.paths.unaligned_reads_paths,
            int(self.args.processes),
            int(self.args.segemehl_accuracy),
            int(self.args.segemehl_evalue))

    def _sam_to_bam(self):
        sam_to_bam_converter = SamToBamConverter()
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self.args.processes) as executor:
            for sam_path, bam_prefix_path in zip(
                    self.paths.read_alignment_result_sam_paths,
                    self.paths.read_alignment_result_bam_prefixes_paths):
                jobs.append(executor.submit(
                    sam_to_bam_converter.sam_to_bam, sam_path, bam_prefix_path))
        # Evaluate thread outcome
        self._check_job_completeness(jobs)

    def _generate_read_alignment_stats(self):
        raw_stat_data_writer = RawStatDataWriter(pretty=True)
        read_files_and_jobs = {}
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self.args.processes) as executor:
            for (read_file, read_alignment_result_bam_path,
                 unaligned_reads_path) in zip(
                     self.read_files,
                     self.paths.read_alignment_result_bam_paths,
                     self.paths.unaligned_reads_paths):
                read_aligner_stats = ReadAlignerStats()
                read_files_and_jobs[read_file]  = executor.submit(
                    read_aligner_stats.count, read_alignment_result_bam_path,
                    unaligned_reads_path)
        # Evaluate thread outcome
        self._check_job_completeness(read_files_and_jobs.values())
        read_files_and_stats = dict(
            [(read_file, job.result()) for read_file, job in
             read_files_and_jobs.items()])
        raw_stat_data_writer.write(
            read_files_and_stats, self.paths.read_aligner_stats_path)

    def _write_alignment_stat_table(self):
        raw_stat_data_reader = RawStatDataReader()
        read_processing_stats = raw_stat_data_reader.read(
            self.paths.read_processing_stats_path)
        alignment_stats = raw_stat_data_reader.read(
            self.paths.read_aligner_stats_path)
        table = []
        table.append(["Lib"] + self.read_files)
        ref_ids = sorted(list(list(alignment_stats.values())[0][
            "countings_per_reference"].keys()))
        table = [
            ["Library read file"] + self.read_files,
            ["No. of input reads"] +
            self._get_read_process_numbers(
                read_processing_stats, "total_no_of_reads"),
            [ "No. of reads - PolyA detected and removed"] +
            self._get_read_process_numbers(
                read_processing_stats, "polya_removed"),
            ["No. of reads - Single 3' A removed"] +
            self._get_read_process_numbers(
                read_processing_stats, "single_a_removed"),
            ["No. of reads - Unmodified"] + self._get_read_process_numbers(
            read_processing_stats, "unmodified"),
            ["No. of reads - Removed as too short"] +
            self._get_read_process_numbers(read_processing_stats, "too_short"),
            ["No. of reads - Long enough and used for alignment"] +
            self._get_read_process_numbers(
                read_processing_stats, "long_enough"),
            ["Total no. of aligned reads"] + [
                round(num) for num in self._total_alignment_stat_numbers(
                alignment_stats, "no_of_aligned_reads")],
            ["Total no. of unaligned reads"] + [str(alignment_stats[read_file][
                "no_of_unaligned_reads"]) for read_file in self.read_files],
            ["Total no. of uniquely aligned reads"] +
            self._total_alignment_stat_numbers(
                alignment_stats, "no_of_uniquely_aligned_reads"),
            ["Total no. of alignments"] + self._total_alignment_stat_numbers(
                alignment_stats, "no_of_alignments"),
            ["Percentage of aligned reads (compared to total input reads)"]  + [
                round(self._calc_percentage(aligned_reads, total_reads), 2)
                for aligned_reads, total_reads in
                zip(self._total_alignment_stat_numbers(
                    alignment_stats, "no_of_aligned_reads"),
                    self._get_read_process_numbers(
                        read_processing_stats, "total_no_of_reads"))],
            ["Percentage of uniquely aligned reads (in relation to all aligned reads)"]  + [
                round(self._calc_percentage(uniquely_aligned_reads, aligned_reads), 2)
                for uniquely_aligned_reads, aligned_reads  in zip(
                        self._total_alignment_stat_numbers(alignment_stats, "no_of_uniquely_aligned_reads"),
                        self._total_alignment_stat_numbers(alignment_stats, "no_of_aligned_reads"))]
        ]
        for ref_id in ref_ids:
            table.append(
                ["%s - No. of aligned reads" % ref_id] +
                self._alignment_number_per_ref_seq(
                    alignment_stats, ref_id, "no_of_aligned_reads"))
            table.append(
                ["%s - No. of uniquely aligned reads" % ref_id] +
                self._alignment_number_per_ref_seq(
                    alignment_stats, ref_id, "no_of_uniquely_aligned_reads"))
            table.append(
                ["%s - No. of alignments" % ref_id] +
                self._alignment_number_per_ref_seq(
                    alignment_stats, ref_id, "no_of_alignments"))
        table_fh = open(self.paths.read_alignment_stats_table_path, "w")
        table_fh.write("\n".join(["\t".join([str(cell) for cell in row]) for row in table]))
        table_fh.close()

    def _calc_percentage(self, mult, div):
        try:
            return(float(mult)/float(div)*100)
        except ZeroDivisionError:
            return(0.0)

    def _alignment_number_per_ref_seq(self, alignment_stats, ref_id, attribute):
        return([alignment_stats[read_file]["countings_per_reference"][
            ref_id][attribute] for read_file in self.read_files])

    def _total_alignment_stat_numbers(self, alignment_stats, attribute):
        return([alignment_stats[read_file]["countings_total"][attribute]
                for read_file in self.read_files])

    def _get_read_process_numbers(
            self, read_processing_stats, attribute):
        return([read_processing_stats[read_file][attribute]
                for read_file in self.read_files])

    def _ref_ids_to_file(self, ref_seq_paths):
        ref_ids_to_file = {}
        fasta_parser = FastaParser()
        for ref_seq_path in ref_seq_paths:
            ref_seq_file = os.path.basename(ref_seq_path)
            ref_seq_id = fasta_parser.header_id(
                fasta_parser.single_entry_file_header(open(ref_seq_path)))
            ref_ids_to_file[ref_seq_id] = ref_seq_file
        return(ref_ids_to_file)

    def create_coverage_files(self):
        """Create coverage files based on the read alignments."""
        read_files = self.paths._get_read_files()
        self.paths.set_read_files_dep_file_lists(read_files)
        raw_stat_data_reader = RawStatDataReader()
        alignment_stats = [
            raw_stat_data_reader.read(
            self.paths.read_aligner_stats_path)]
        read_files_and_no_of_aligned_reads = dict([
            (read_file,
             round(attributes["countings_total"]["no_of_aligned_reads"]))
             for read_file, attributes in alignment_stats[0].items()])
        min_read_alignment_counting = min(
            read_files_and_no_of_aligned_reads.values())
        # Run the generation of coverage in parallel
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self.args.processes) as executor:
            for read_file, bam_path in zip(
                read_files, self.paths.read_alignment_result_bam_paths):
                jobs.append(executor.submit(
                        self._create_coverage_files_for_lib,
                        read_file, bam_path,
                        read_files_and_no_of_aligned_reads,
                        min_read_alignment_counting))
        # Evaluate thread outcome
        self._check_job_completeness(jobs)

    def _create_coverage_files_for_lib(
        self, read_file, bam_path, read_alignment_stats,
        min_read_alignment_counting):
        coverage_creator = CoverageCreator()
        read_count_splitting = True
        if self.args.skip_read_count_splitting is True:
            read_count_splitting = False
        coverage_creator.init_coverage_lists(bam_path)
        coverage_creator.count_coverage(
            bam_path, read_count_splitting=read_count_splitting,
            uniqueley_aligned_only=self.args.unique_only,
            first_base_only=self.args.first_base_only)
        # Raw countings
        coverage_creator.write_to_files(
            "%s/%s" % (self.paths.coverage_folder, read_file),
            read_file)
        total_number_of_aligned_reads = read_alignment_stats[read_file]
        # Read normalized countings - multiplied by min read counting
        factor = (min_read_alignment_counting / total_number_of_aligned_reads)
        coverage_creator.write_to_files(
            "%s/%s-div_by_%.1f_multi_by_%.1f" % (
                self.paths.coverage_folder_norm_reads, read_file,
                total_number_of_aligned_reads, min_read_alignment_counting),
            read_file, factor=factor)
        # Read normalized countings - multiplied by 1M
        factor = (1000000 / total_number_of_aligned_reads)
        coverage_creator.write_to_files(
            "%s/%s-div_by_%.1f_multi_by_1M" % (
                self.paths.coverage_folder_norm_reads_mil, read_file,
                total_number_of_aligned_reads), read_file, factor=factor)

    def _check_job_completeness(self, jobs):
        """Check the completness of each job in a list"""
        for job in concurrent.futures.as_completed(jobs):
            if job.exception():
                raise(job.exception())

    def quantify_gene_wise(self):
        norm_by_alignment_freq = True
        norm_by_overlap_freq = True
        if self.args.skip_norm_by_alignment_freq:
            norm_by_alignment_freq = False
        if self.args.skip_norm_by_overlap_freq:
            norm_by_overlap_freq = False
        read_files = self.paths._get_read_files()
        annotation_files = self.paths._get_annotation_files()
        self.paths.set_annotation_paths(annotation_files)
        self.paths.set_read_files_dep_file_lists(read_files)
        for read_file, read_alignment_path in zip(
                read_files, self.paths.read_alignment_result_bam_paths):
            gene_wise_quantification = GeneWiseQuantification(
                min_overlap=self.args.min_overlap,
                norm_by_alignment_freq=norm_by_alignment_freq,
                norm_by_overlap_freq=norm_by_overlap_freq)
            gene_wise_quantification.calc_overlaps_per_alignment(
                read_alignment_path, self.paths.annotation_paths)
            for  annotation_file, annotation_path in zip(
                    annotation_files, self.paths.annotation_paths):
                gene_wise_quantification.quantify(
                    read_alignment_path, annotation_path,
                    self.paths.gene_quanti_path(
                        read_file, annotation_file))
        self._gene_quanti_create_overview(
            annotation_files, self.paths.annotation_paths, read_files)

    def _gene_quanti_create_overview(
            self, annotation_files, annotation_paths, read_files):
        gene_wise_overview = GeneWiseOverview()
        path_and_name_combos = {}
        for annotation_file, annotation_path in zip(
                annotation_files, annotation_paths):
            path_and_name_combos[annotation_path] = []
            for read_file in read_files:
                path_and_name_combos[annotation_path].append(
                    [read_file, self.paths.gene_quanti_path(
                        read_file, annotation_file)])
        gene_wise_overview.create_overview(
            path_and_name_combos, read_files,
            self.paths.gene_wise_quanti_combined_path)
