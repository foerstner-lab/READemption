import concurrent.futures
import os
import sys
import json
sys.path.append(".")
from libs.coveragecalculator import CoverageCalculator
from libs.deseq import DESeqRunner
from libs.fasta import FastaParser
from libs.genewisequanti import GeneWiseQuantification, GeneWiseOverview
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

    """Manage the actions of the subcommands.

    The Controller take care of providing the argumentes like path
    names and the parallel processing of tasks.

    """

    def __init__(self, args):
        """Create an instance."""
        self._args = args
        self._paths = Paths(args.project_path)
        self._read_files = None
        self._ref_seq_files = None
        self._lib_names = None

    def create_project(self, version):
        """Create a new project."""
        project_creator = ProjectCreator()
        project_creator.create_root_folder(self._args.project_path)
        project_creator.create_subfolders(self._paths.required_folders())
        project_creator.create_version_file(self._paths.version_path, version)
        sys.stdout.write("Created folder \"%s\" and required subfolders.\n" % (
                self._args.project_path))
        sys.stdout.write("Please copy read files into folder \"%s\" and "
                         "reference sequences files into folder \"%s\".\n" % (
                self._paths.read_fasta_folder, self._paths.ref_seq_folder))

    def align_reads(self):
        """Perform the alignment of the reads."""
        self._test_folder_existance(
            self._paths.required_read_alignment_folders())
        self._read_files = self._paths.get_read_files()
        self._lib_names = self._paths.get_lib_names()
        self._ref_seq_files = self._paths.get_ref_seq_files()
        self._test_align_file_existance()
        self._paths.set_read_files_dep_file_lists(
            self._read_files, self._lib_names)
        self._paths.set_ref_seq_paths(self._ref_seq_files)
        self._prepare_reads()
        self._align_reads()
        self._sam_to_bam()
        self._generate_read_alignment_stats()
        self._write_alignment_stat_table()

    def _test_align_file_existance(self):
        """Test if the input file for the the align subcommand exist."""
        if len(self._read_files) == 0:
            self._write_err_msg_and_quit("Error! No read libraries given!\n")
        if len(self._ref_seq_files ) == 0:
            self._write_err_msg_and_quit("Error! No reference sequence files given!\n")
        
    def _test_folder_existance(self, task_specific_folders):
        """Test the existance of required folders."""
        for folder in (
            self._paths.required_base_folders() + task_specific_folders):
            if not os.path.exists(folder):
                self._write_err_msg_and_quit(
                    "Error! Folder '%s' does not exist! Is the given project "
                    "folder name correct?\n" % folder)

    def _file_needs_to_be_created(self, file_path, quiet=False):
        """Test if a file exists of need to be created."""
        if self._args.force is True:
            return True
        if os.path.exists(file_path):
            if quiet is False:
                sys.stderr.write(
                    "File %s exists. Skipping its generation.\n" % file_path)
            return False
        return True

    def _prepare_reads(self):
        """Manage the prepartion of reads before the actual mappings."""
        raw_stat_data_writer = RawStatDataWriter(pretty=True)
        read_files_and_jobs = {}
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes) as executor:
            for lib_name, read_path, processed_read_path in zip(
                    self._lib_names, self._paths.read_paths, 
                    self._paths.processed_read_paths):
                if self._file_needs_to_be_created(processed_read_path) is False:
                    continue
                read_processor = ReadProcessor(
                    poly_a_clipping=self._args.poly_a_clipping,
                    min_read_length=self._args.min_read_length)
                read_files_and_jobs[lib_name]  = executor.submit(
                    read_processor.process, read_path, processed_read_path)
        # Evaluate thread outcome
        self._check_job_completeness(read_files_and_jobs.values())
        if self._file_needs_to_be_created(
            self._paths.read_processing_stats_path) is False:
            return
        # Create a dict of the read file names and the processing
        # counting results
        read_files_and_stats = dict(
            [(lib_name, job.result()) for lib_name, job in
             read_files_and_jobs.items()])
        raw_stat_data_writer.write(
            read_files_and_stats, self._paths.read_processing_stats_path)

    def _align_reads(self):
        """Manage the alignemnt of reads."""
        read_aligner = ReadAligner(self._args.segemehl_bin, self._args.progress)
        if self._file_needs_to_be_created(self._paths.index_path) is True:
            read_aligner.build_index(
                self._paths.ref_seq_paths, self._paths.index_path)
        for read_path, output_path, nomatch_path, bam_path in zip(
            self._paths.processed_read_paths, 
            self._paths.read_alignment_result_sam_paths, 
            self._paths.unaligned_reads_paths, 
            self._paths.read_alignment_result_bam_paths):
            if self._file_needs_to_be_created(output_path) is False:
                continue
            elif self._file_needs_to_be_created(bam_path) is False:
                continue
            read_aligner.run_alignment(
                read_path, self._paths.index_path, self._paths.ref_seq_paths, 
                output_path, nomatch_path, int(self._args.processes),
                int(self._args.segemehl_accuracy), 
                float(self._args.segemehl_evalue), self._args.split)

    def _sam_to_bam(self):
        """Manage the conversion of mapped read from SAM to BAM format."""
        sam_to_bam_converter = SamToBamConverter()
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes) as executor:
            for sam_path, bam_prefix_path, bam_path in zip(
                    self._paths.read_alignment_result_sam_paths,
                    self._paths.read_alignment_result_bam_prefixes_paths,
                    self._paths.read_alignment_result_bam_paths):
                if self._file_needs_to_be_created(bam_path) is False:
                    continue
                jobs.append(executor.submit(
                    sam_to_bam_converter.sam_to_bam, sam_path, bam_prefix_path))
        # Evaluate thread outcome
        self._check_job_completeness(jobs)

    def _generate_read_alignment_stats(self):
        """Manage the generation of alingment statistics."""
        raw_stat_data_writer = RawStatDataWriter(pretty=True)
        read_files_and_jobs = {}
        if self._file_needs_to_be_created(
            self._paths.read_aligner_stats_path) is False:
            return
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes) as executor:
            for (lib_name, read_alignment_result_bam_path,
                 unaligned_reads_path) in zip(
                     self._lib_names, 
                     self._paths.read_alignment_result_bam_paths,
                     self._paths.unaligned_reads_paths):
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
            read_files_and_stats, self._paths.read_aligner_stats_path)

    def _write_alignment_stat_table(self):
        """Manage the creation of the mapping statistic output table."""
        raw_stat_data_reader = RawStatDataReader()
        read_processing_stats = raw_stat_data_reader.read(
            self._paths.read_processing_stats_path)
        alignment_stats = raw_stat_data_reader.read(
            self._paths.read_aligner_stats_path)
        read_aligner_stats_table = ReadAlignerStatsTable(
            read_processing_stats, alignment_stats, 
            self._lib_names, self._paths.read_alignment_stats_table_path)
        read_aligner_stats_table.write()

    def _ref_ids_to_file(self, ref_seq_paths):
        """Translate the reference ID to file paths."""
        ref_ids_to_file = {}
        fasta_parser = FastaParser()
        for ref_seq_path in ref_seq_paths:
            ref_seq_file = os.path.basename(ref_seq_path)
            with open(ref_seq_path) as ref_seq_fh:
                ref_seq_id = fasta_parser.header_id(
                    fasta_parser.single_entry_file_header(ref_seq_fh))
                ref_ids_to_file[ref_seq_id] = ref_seq_file
        return ref_ids_to_file

    def create_coverage_files(self):
        """Create coverage files based on the read alignments.

        The coverages are calculated per replicon and the results are
        written to the output file. This might be slower but if all
        coverages are detmined at once the data structure will become
        too large when working with large reference sequences.

        """
        self._test_folder_existance(self._paths.required_coverage_folders())
        lib_names = self._paths.get_lib_names()
        self._paths.set_read_files_dep_file_lists(
            self._paths.get_read_files(), lib_names)
        raw_stat_data_reader = RawStatDataReader()
        alignment_stats = [raw_stat_data_reader.read(
                self._paths.read_aligner_stats_path)]
        # Get number of aligned of number of uniquely aligned reads
        if self._args.unique_only is False:
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
                lib_names, self._paths.read_alignment_result_bam_paths):
                no_of_aligned_reads = float(
                    read_files_aligned_read_freq[lib_name])
                jobs.append(executor.submit(
                        self._create_coverage_files_for_lib,
                        lib_name, bam_path, no_of_aligned_reads,
                        min_no_of_aligned_reads))
        # Evaluate thread outcome
        self._check_job_completeness(jobs)

    def _all_coverage_file_exist(
        self, lib_name, strands, no_of_aligned_reads, min_no_of_aligned_reads):
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
        if any([self._file_needs_to_be_created(file, quiet=True) 
                    for file in files]) is False:
            sys.stderr.write(
                "The files %s exists. Skipping their generation.\n" % 
                ", " .join(files))
            return True
        return False

    def _create_coverage_files_for_lib(
        self, lib_name, bam_path, no_of_aligned_reads, min_no_of_aligned_reads):
        """Perform the coverage calculation for a given library."""
        strands = ["forward", "reverse"]
        if self._all_coverage_file_exist(
            lib_name, strands, no_of_aligned_reads, min_no_of_aligned_reads):
            return
        read_count_splitting = True
        if self._args.skip_read_count_splitting is True:
            read_count_splitting = False
        coverage_calculator = CoverageCalculator(
            read_count_splitting=read_count_splitting,
            uniqueley_aligned_only=self._args.unique_only,
            first_base_only=self._args.first_base_only)
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
                                div=no_of_aligned_reads), "w")))
                                                for strand in strands])
        coverage_writers_tnoar_mil_norm = dict([(
                    strand, WiggleWriter(
                        "%s_%s" % (lib_name, strand),
                        open(self._paths.wiggle_file_tnoar_norm_mil_path(
                                lib_name, strand, multi=1000000,
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
        """Manage the counting of alinged read per gene."""
        self._test_folder_existance(
            self._paths.required_gene_quanti_folders())
        norm_by_alignment_freq = True
        norm_by_overlap_freq = True
        if self._args.skip_norm_by_alignment_freq:
            norm_by_alignment_freq = False
        if self._args.skip_norm_by_overlap_freq:
            norm_by_overlap_freq = False
        lib_names = self._paths.get_lib_names()
        annotation_files = self._paths.get_annotation_files()
        self._paths.set_annotation_paths(annotation_files)
        self._paths.set_read_files_dep_file_lists(
            self._paths.get_read_files(), lib_names)
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes) as executor:
            for lib_name, read_alignment_path in zip(
                lib_names, self._paths.read_alignment_result_bam_paths):
                jobs.append(executor.submit(
                        self._quantify_gene_wise, lib_name,
                        read_alignment_path, norm_by_alignment_freq,
                        norm_by_overlap_freq, annotation_files))
        # Evaluate thread outcome
        self._check_job_completeness(jobs)
        self._gene_quanti_create_overview(
            annotation_files, self._paths.annotation_paths, lib_names)

    def _quantify_gene_wise(self, lib_name, read_alignment_path, 
                    norm_by_alignment_freq,  norm_by_overlap_freq, 
                    annotation_files):
        """Perform the gene wise quantification for a given library."""
        gene_quanti_paths = [
            self._paths.gene_quanti_path(lib_name, annotation_file)
            for annotation_file in annotation_files]
        # Check if all output files for this library exist - if so
        # skip their creation
        if any([self._file_needs_to_be_created(gene_quanti_path, quiet=True)
                for gene_quanti_path in gene_quanti_paths]) is False:
            sys.stderr.write(
                "The file(s) %s exist(s). Skipping their/its generation.\n" % 
                ", " .join(gene_quanti_paths))
            return
        gene_wise_quantification = GeneWiseQuantification(
            min_overlap=self._args.min_overlap,
            norm_by_alignment_freq=norm_by_alignment_freq,
            norm_by_overlap_freq=norm_by_overlap_freq,
            allowed_features_str=self._args.allowed_features,
            skip_antisense=self._args.skip_antisense,
            unique_only=self._args.unique_only)
        gene_wise_quantification.calc_overlaps_per_alignment(
            read_alignment_path, self._paths.annotation_paths)
        for annotation_file, annotation_path in zip(
            annotation_files, self._paths.annotation_paths):
            gene_wise_quantification.quantify(
                read_alignment_path, annotation_path,
                self._paths.gene_quanti_path(
                    lib_name, annotation_file), self._args.pseudocounts)

    def _gene_quanti_create_overview(
            self, annotation_files, annotation_paths, lib_names):
        """Create na overview table of all gene quantification for all libs."""
        gene_wise_overview = GeneWiseOverview(
            allowed_features_str=self._args.allowed_features,
            skip_antisense=self._args.skip_antisense)
        path_and_name_combos = {}
        for annotation_file, annotation_path in zip(
                annotation_files, annotation_paths):
            path_and_name_combos[annotation_path] = []
            for read_file in lib_names:
                path_and_name_combos[annotation_path].append(
                    [read_file, self._paths.gene_quanti_path(
                        read_file, annotation_file)])
        if self._file_needs_to_be_created(
            self._paths.gene_wise_quanti_combined_path) is True:
            gene_wise_overview.create_overview_raw_countings(
                path_and_name_combos, lib_names,
                self._paths.gene_wise_quanti_combined_path)
        if self._file_needs_to_be_created(
            self._paths.gene_wise_quanti_combined_rpkm_path) is True:
            gene_wise_overview.create_overview_rpkm(
                path_and_name_combos, lib_names,
                self._paths.gene_wise_quanti_combined_rpkm_path,
                self._libs_and_total_num_of_aligned_reads())
        if self._file_needs_to_be_created(
            self._paths.gene_wise_quanti_combined_tnoar_path) is True:
            gene_wise_overview.create_overview_norm_by_tnoar(
                path_and_name_combos, lib_names,
                self._paths.gene_wise_quanti_combined_tnoar_path,
                self._libs_and_total_num_of_aligned_reads())

    def _libs_and_total_num_of_aligned_reads(self):
        """Read the total number of reads per library."""
        with open(self._paths.read_aligner_stats_path) as read_aligner_stats_fh:
            read_aligner_stats = json.loads(read_aligner_stats_fh.read())
        return dict([(lib, values["stats_total"]["no_of_aligned_reads"])
                     for lib, values in read_aligner_stats.items()])

    def _libs_and_total_num_of_uniquely_aligned_reads(self):
        """Read the total number of reads per library."""
        with open(self._paths.read_aligner_stats_path) as read_aligner_stats_fh:
            read_aligner_stats = json.loads(read_aligner_stats_fh.read())
        return dict([(lib, values["stats_total"]["no_of_uniquely_aligned_reads"])
                     for lib, values in read_aligner_stats.items()])

    def compare_with_deseq(self):
        """Manage the pairwise expression comparison with DESeq."""
        self._test_folder_existance(
            self._paths.required_deseq_folders())
        arg_libs = [self._paths._clean_file_name(lib) for lib in 
                self._args.libs.split(",")]
        conditions = self._args.conditions.split(",")
        self._check_deseq_args(arg_libs, conditions)
        deseq_runner = DESeqRunner(
            arg_libs, conditions, self._paths.deseq_raw_folder,
            self._paths.deseq_extended_folder, self._paths.deseq_script_path,
            self._paths.gene_wise_quanti_combined_path, 
            self._paths.deseq_tmp_session_info_script,
            self._paths.deseq_session_info,
            no_replicates=self._args.no_replicates)
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
                "(= %s elements)\nand \n%s (= %s elements).\n" % (
                    self._args.libs, len(arg_libs), self._args.conditions,
                    len(conditions)))
        libs = self._paths.get_lib_names()
        if len(libs) != len(arg_libs):
            self._write_err_msg_and_quit(
                "The number of read libraries is lower or higher than "
                "expected. The following read libs are available: %s\nThe "
                "following read list string is suggested: \"%s\"\n" % (
                    ", ".join(read_files), ",".join(libs)))
        for lib in libs:
            if lib not in arg_libs:
                self._write_err_msg_and_quit(
                    "The library \"%s\" is not present in your list of "
                    "libraries. Please add it.\n" % (lib))

    def _write_err_msg_and_quit(self, msg):
        """Write error message and close the program gracefully."""
        sys.stderr.write(msg)
        sys.exit(1)

    def viz_align(self):
        """Generate plots based on the read processing and mapping"""
        from libs.vizalign import AlignViz
        align_viz = AlignViz(
            self._paths.get_lib_names(),
            self._paths.read_processing_stats_path,
            self._paths.read_aligner_stats_path)
        align_viz.read_stat_files()
        align_viz.plot_input_read_length(
            self._paths.viz_align_input_read_length_plot_path)
        align_viz.plot_processed_read_length(
            self._paths.viz_align_processed_reads_length_plot_path)

    def viz_gene_quanti(self):
        """Generate plots based on the gene-wise read countings"""
        from libs.vizgenequanti import GeneQuantiViz
        gene_quanti_viz = GeneQuantiViz(
            self._paths.gene_wise_quanti_combined_path, 
            self._paths.get_lib_names())
        gene_quanti_viz.parse_input_table()
        gene_quanti_viz.plot_correlations(
            self._paths.viz_gene_quanti_scatter_plot_path)
        gene_quanti_viz.plot_annotation_class_quantification(
            self._paths.viz_gene_quanti_rna_classes_plot_path)

    def viz_deseq(self):
        """Generate plots based on the DESeq analysis"""
        from libs.vizdeseq import DESeqViz
        deseq_path_template = (
            self._paths.deseq_raw_folder + "/deseq_comp_%s_vs_%s.csv")
        deseq_viz = DESeqViz(self._paths.deseq_script_path, deseq_path_template)
        deseq_viz.create_scatter_plots(
            self._paths.viz_deseq_scatter_plot_path)
        deseq_viz.create_volcano_plots(
            self._paths.viz_deseq_volcano_plot_path,
            self._paths.viz_deseq_volcano_plot_adj_path)

