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
        self.report_folder = "%s/reports_and_stats" % (self.output_folder)
        self.raw_stat_data_folder = "%s/stats_data_json" % (self.output_folder)
        self._set_input_folder_names()
        self._set_read_alignment_folder_names()
        self._set_coverage_folder_names()
        self._set_gene_quanti_folder_names()

    def _set_input_folder_names(self):
        self.read_fasta_folder = "%s/reads" % self.input_folder
        self.ref_seq_folder = "%s/reference_sequences" % self.input_folder
        self.annotation_folder = "%s/annotation_files" % self.input_folder

    def _set_read_alignment_folder_names(self):
        self.read_alignment_index_folder = "%s/read_alignments-index" % (
            self.output_folder)
        self.read_alignments_folder = "%s/read_alignments-alignments" % (
            self.output_folder)
        self.processed_reads_folder = "%s/read_alignments-processed_reads" % (
            self.output_folder)
        self.unaligned_reads_folder = "%s/read_alignments-unaligned_reads" % (
            self.output_folder)

    def _set_coverage_folder_names(self):
        base_name = "%s/coverages" % self.output_folder
        self.coverage_raw_folder = "%s-raw" % base_name
        self.coverage_tnoar_min_norm_folder = "%s-tnoar_min_normalized" % (
            base_name)
        self.coverage_tnoar_mil_norm_folder = "%s-tnoar_mil_normalized" % (
            base_name)

    def _set_gene_quanti_folder_names(self):
        self.gene_quanti_folder = (
            "%s/gene_wise_quantifications" % self.output_folder)
        self.gene_wise_quanti_combined_path = (
            "%s/Gene_wise_quantifications_combined.csv" % self.output_folder)

    def _set_static_files(self):
        """Set name of common files."""
        self.read_processing_stats_path = "%s/read_processing.json" % (
            self.raw_stat_data_folder)
        self.read_aligner_stats_path = "%s/read_alignment.json" % (
            self.raw_stat_data_folder)
        self.read_file_stats = "%s/input_read_stats.txt" % (
            self.report_folder)
        self.ref_seq_file_stats = "%s/reference_sequences_stats.txt" % (
            self.report_folder)
        self.annotation_file_stats = "%s/annotation_file_stats.txt" % (
            self.report_folder)
        self.read_alignment_stats_table_path = "%s/read_alignment_stats.csv" % (
            self.report_folder)
        self.index_path = "%s/index.idx" % self.read_alignment_index_folder
        self.deseq_file = "%s/deseq.R" % self.output_folder

    def _get_sorted_folder_content(self, folder):
        """Return the sorted file list of a folder"""
        return(list(filter(lambda file: not file.endswith("~"),
                      sorted(os.listdir(folder)))))

    def get_read_files(self):
        """Read the names of the read files."""
        return(self._get_sorted_folder_content(self.read_fasta_folder))

    def get_ref_seq_files(self):
        """Read the names of reference sequence files."""
        return(self._get_sorted_folder_content(self.ref_seq_folder))

    def get_annotation_files(self):
        """Read the names of annotation files."""
        return(self._get_sorted_folder_content(self.annotation_folder))

    def required_folders(self):
        return(self._required_base_folders() +
               self._required_input_folders() +
               self._required_read_alignment_folders() +
               self._required_coverage_folders() +
               self._required_gene_quanti_folder())

    def _required_base_folders(self):
        return([self.input_folder, self.output_folder, self.report_folder,
                self.raw_stat_data_folder])

    def _required_input_folders(self):
        return([self.read_fasta_folder, self.ref_seq_folder,
                self.annotation_folder])

    def _required_read_alignment_folders(self):
        return([self.read_alignments_folder, self.processed_reads_folder,
                self.unaligned_reads_folder, self.read_alignment_index_folder])

    def _required_coverage_folders(self):
        return([self.coverage_raw_folder, self.coverage_tnoar_min_norm_folder,
                self.coverage_tnoar_mil_norm_folder])

    def _required_gene_quanti_folder(self):
        return([self.gene_quanti_folder])

    def set_read_files_dep_file_lists(self, read_files):
        self.read_paths = self._path_list(self.read_fasta_folder, read_files)
        self.processed_read_paths = self._path_list(
            self.processed_reads_folder, read_files, appendix="_processed.fa")
        self.read_alignment_result_sam_paths = self._path_list(
            self.read_alignments_folder, read_files, appendix="_alignments.sam")
        # samtool appends ".bam" so only the prefix is required
        self.read_alignment_result_bam_prefixes_paths = self._path_list(
            self.read_alignments_folder, read_files, appendix="_alignments")
        self.read_alignment_result_bam_paths = self._path_list(
            self.read_alignments_folder, read_files, appendix="_alignments.bam")
        self.unaligned_reads_paths = self._path_list(
            self.unaligned_reads_folder, read_files, appendix="_unaligned.fa")

    def set_ref_seq_paths(self, ref_seq_files):
        self.ref_seq_paths = self._path_list(self.ref_seq_folder, ref_seq_files)

    def set_annotation_paths(self, annotation_files):
        self.annotation_paths = self._path_list(
            self.annotation_folder, annotation_files)

    def _path_list(self, folder, files, appendix=""):
        return(["%s/%s%s" % (folder, file, appendix) for file in files])

    def gene_quanti_path(self, read_file, annotation_file):
        return("%s/%s_to_%s.csv" % (
            self.gene_quanti_folder, read_file, annotation_file))

    def wiggle_file_raw_path(self, read_file, strand, multi=None, div=None):
        return(self._wiggle_file_path(
            self.coverage_raw_folder, read_file, strand, multi=None, div=None))

    def wiggle_file_tnoar_norm_min_path(
            self, read_file, strand, multi=None, div=None):
        return(self._wiggle_file_path(
            self.coverage_tnoar_min_norm_folder, read_file, strand, multi, div))

    def wiggle_file_tnoar_norm_mil_path(
            self, read_file, strand, multi=None, div=None):
        return(self._wiggle_file_path(
            self.coverage_tnoar_mil_norm_folder, read_file, strand, multi, div))

    def _wiggle_file_path(
            self, folder, read_file, strand, multi=None, div=None):
        path = "%s/%s" % (folder, read_file)
        if not div is None:
            path += "_div_by_%.1f" % (div)
        if not multi is None:
            path += "_multi_by_%.1f" % (multi)
        path += "_%s.wig" % strand
        return(path)
