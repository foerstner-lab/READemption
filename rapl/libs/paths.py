import os

class Paths(object):

    def __init__(self, base_path="."):
        self.base_path = base_path
        self._set_folder_names()
        self._set_static_files()

    def _set_folder_names(self):
        """Set the name of folders used in a project."""
        self.input_folder = "%s/input" % (self.base_path)
        self.output_folder = "%s/output" % (self.base_path)
        self._set_input_folder_names()
        self._set_read_alignment_folder_names()
        self._set_coverage_folder_names()
        self._set_gene_quanti_folder_names()
        self._set_deseq_folder_names()
        self._set_viz_gene_quanti_folder_names()
        self._set_viz_deseq_folder_names()

    def _set_input_folder_names(self):
        self.read_fasta_folder = "%s/reads" % self.input_folder
        self.ref_seq_folder = "%s/reference_sequences" % self.input_folder
        self.annotation_folder = "%s/annotation_files" % self.input_folder

    def _set_read_alignment_folder_names(self):
        self.align_base_folder = "%s/align" % self.output_folder
        self.read_alignment_index_folder = "%s/index" % (
            self.align_base_folder)
        self.read_alignments_folder = "%s/alignments" % (
            self.align_base_folder)
        self.processed_reads_folder = "%s/processed_reads" % (
            self.align_base_folder)
        self.unaligned_reads_folder = "%s/unaligned_reads" % (
            self.align_base_folder)
        self.align_report_folder = "%s/reports_and_stats" % (
            self.align_base_folder)
        self.raw_stat_data_folder = "%s/stats_data_json" % (
            self.align_report_folder)
        self.raw_stat_data_folder = "%s/stats_data_json" % (
            self.align_report_folder)

    def _set_coverage_folder_names(self):
        self.coverage_base_folder = "%s/coverage" % self.output_folder
        self.coverage_raw_folder = "%s/coverage-raw" % self.coverage_base_folder
        self.coverage_tnoar_min_norm_folder = "%s/coverage-tnoar_min_normalized" % (
            self.coverage_base_folder)
        self.coverage_tnoar_mil_norm_folder = "%s/coverage-tnoar_mil_normalized" % (
            self.coverage_base_folder)

    def _set_gene_quanti_folder_names(self):
        self.gene_quanti_base_folder = (
            "%s/gene_quanti" % self.output_folder)
        self.gene_quanti_per_lib_folder = "%s/gene_quanti_per_lib" % (
            self.gene_quanti_base_folder)
        self.gene_quanti_combined_folder = "%s/gene_quanti_combined" % (
            self.gene_quanti_base_folder)
        self.gene_wise_quanti_combined_path = (
            "%s/gene_wise_quantifications_combined.csv" % 
            self.gene_quanti_combined_folder)
        self.gene_wise_quanti_combined_rpkm_path = (
            "%s/gene_wise_quantifications_combined_rpkm.csv" % 
            self.gene_quanti_combined_folder)
        self.gene_wise_quanti_combined_tnoar_path = (
            "%s/gene_wise_quantifications_combined_tnoar.csv" % 
            self.gene_quanti_combined_folder)

    def _set_deseq_folder_names(self):
        self.deseq_base_folder = ("%s/deseq" % self.output_folder)
        self.deseq_raw_folder = ("%s/deseq_raw" % self.deseq_base_folder)
        self.deseq_extended_folder = (
            "%s/deseq_with_annotations" % self.deseq_base_folder)

    def _set_viz_gene_quanti_folder_names(self):
        self.viz_gene_quanti_base_folder = (
            "%s/viz_gene_quanti" % self.output_folder)
        self.viz_gene_quanti_scatter_plot_path = (
            "%s/expression_scatter_plots.pdf" % 
            self.viz_gene_quanti_base_folder)

    def _set_viz_deseq_folder_names(self):
        self.viz_deseq_base_folder = (
            "%s/viz_deseq" % self.output_folder)
        self.viz_deseq_volcano_plot_path = (
            "%s/volcano_plots_log2_fold_change_vs_p-value.pdf" % 
            self.viz_deseq_base_folder)
        self.viz_deseq_volcano_plot_adj_path = (
            "%s/volcano_plots_log2_fold_change_vs_adjusted_p-value.pdf" %
            self.viz_deseq_base_folder)

    def _set_static_files(self):
        """Set name of common files."""
        self.read_processing_stats_path = "%s/read_processing.json" % (
            self.raw_stat_data_folder)
        self.read_aligner_stats_path = "%s/read_alignment.json" % (
            self.raw_stat_data_folder)
        self.read_file_stats = "%s/input_read_stats.txt" % (
            self.align_report_folder)
        self.ref_seq_file_stats = "%s/reference_sequences_stats.txt" % (
            self.align_report_folder)
        #self.annotation_file_stats = "%s/annotation_file_stats.txt" % (
        #    self.report_folder)
        self.read_alignment_stats_table_path = "%s/read_alignment_stats.csv" % (
            self.align_report_folder)
        self.index_path = "%s/index.idx" % self.read_alignment_index_folder
        self.deseq_script_path = "%s/deseq.R" % self.deseq_raw_folder
        self.deseq_tmp_session_info_script = "%s/tmp.R" % self.deseq_raw_folder
        self.deseq_session_info = "%s/R_session_info.txt" % (
            self.deseq_raw_folder)
        self.version_path = "%s/used_rapl_version.txt" % (self.align_report_folder)

    def _get_sorted_folder_content(self, folder):
        """Return the sorted file list of a folder"""
        return list(filter(lambda file: 
                           not (file.endswith("~") or 
                                os.path.basename(file).startswith(".")),
                      sorted(os.listdir(folder))))

    def get_read_files(self):
        """Read the names of the read files."""
        return self._get_sorted_folder_content(self.read_fasta_folder)

    def get_lib_names(self):
        return [self._clean_file_name(file_name) 
                for file_name in self.get_read_files()]

    def _clean_file_name(self, file_name):
        for suffix in ["bz2", "BZ", "gz", "GZ", "fa", "fasta", "FA", "FASTA"]:
            if file_name.endswith(suffix):
                suffix = "." + suffix
                file_name = file_name[:-len(suffix)]
        return file_name

    def get_ref_seq_files(self):
        """Read the names of reference sequence files."""
        return self._get_sorted_folder_content(self.ref_seq_folder)

    def get_annotation_files(self):
        """Read the names of annotation files."""
        return self._get_sorted_folder_content(self.annotation_folder)

    def required_folders(self):
        return (self.required_base_folders() +
                self.required_input_folders() +
                self.required_read_alignment_folders() +
                self.required_coverage_folders() +
                self.required_gene_quanti_folders() +
                self.required_deseq_folders() +
                self.required_viz_gene_quanti_folders() + 
                self.required_viz_deseq_folders())

    def required_base_folders(self):
        return [self.input_folder, self.output_folder]

    def required_input_folders(self):
        return [self.read_fasta_folder, self.ref_seq_folder,
                self.annotation_folder]

    def required_read_alignment_folders(self):
        return [self.align_base_folder, self.read_alignments_folder, 
                self.processed_reads_folder, self.unaligned_reads_folder, 
                self.read_alignment_index_folder, self.align_report_folder,
                self.raw_stat_data_folder]

    def required_coverage_folders(self):
        return [self.coverage_base_folder, self.coverage_raw_folder, 
                self.coverage_tnoar_min_norm_folder, 
                self.coverage_tnoar_mil_norm_folder]

    def required_gene_quanti_folders(self):
        return [self.gene_quanti_base_folder, self.gene_quanti_per_lib_folder, 
                self.gene_quanti_combined_folder]

    def required_deseq_folders(self):
        return [self.deseq_base_folder, self.deseq_raw_folder, 
                self.deseq_extended_folder]

    def required_viz_gene_quanti_folders(self):
        return [self.viz_gene_quanti_base_folder]

    def required_viz_deseq_folders(self):
        return [self.viz_deseq_base_folder]

    def set_read_files_dep_file_lists(self, read_files, lib_names):
        self.read_paths = self._path_list(self.read_fasta_folder, read_files)
        self.processed_read_paths = self._path_list(
            self.processed_reads_folder, lib_names,
            appendix="_processed.fa.gz")
        self.read_alignment_result_sam_paths = self._path_list(
            self.read_alignments_folder, lib_names,
            appendix="_alignments.sam")
        # samtool appends ".bam" so only the prefix is required
        self.read_alignment_result_bam_prefixes_paths = self._path_list(
            self.read_alignments_folder, 
            lib_names, appendix="_alignments")
        self.read_alignment_result_bam_paths = self._path_list(
            self.read_alignments_folder, lib_names, 
            appendix="_alignments.bam")
        self.unaligned_reads_paths = self._path_list(
            self.unaligned_reads_folder, lib_names, 
            appendix="_unaligned.fa")

    def set_ref_seq_paths(self, ref_seq_files):
        self.ref_seq_paths = self._path_list(self.ref_seq_folder, ref_seq_files)

    def set_annotation_paths(self, annotation_files):
        self.annotation_paths = self._path_list(
            self.annotation_folder, annotation_files)

    def _path_list(self, folder, files, appendix=""):
        return ["%s/%s%s" % (folder, file, appendix) for file in files]

    def gene_quanti_path(self, read_file, annotation_file):
        return "%s/%s_to_%s.csv" % (
            self.gene_quanti_per_lib_folder, read_file, annotation_file)

    def wiggle_file_raw_path(self, read_file, strand, multi=None, div=None):
        return self._wiggle_file_path(
            self.coverage_raw_folder, read_file, strand, multi=None, div=None)

    def wiggle_file_tnoar_norm_min_path(
            self, read_file, strand, multi=None, div=None):
        return self._wiggle_file_path(
            self.coverage_tnoar_min_norm_folder, read_file, strand, multi, div)

    def wiggle_file_tnoar_norm_mil_path(
            self, read_file, strand, multi=None, div=None):
        return self._wiggle_file_path(
            self.coverage_tnoar_mil_norm_folder, read_file, strand, multi, div)

    def _wiggle_file_path(
            self, folder, read_file, strand, multi=None, div=None):
        path = "%s/%s" % (folder, read_file)
        if not div is None:
            path += "_div_by_%.1f" % (div)
        if not multi is None:
            path += "_multi_by_%.1f" % (multi)
        path += "_%s.wig" % strand
        return path
