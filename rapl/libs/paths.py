import os

class Paths(object):

    def __init__(self, base_path="."):
        self.base_path = base_path
        self._set_folder_names()
        self._set_static_file_names()

    def _set_folder_names(self):
        """Set the name of folders used in a project."""
        self.input_folder = "%s/input" % (self.base_path)
        self.output_folder = "%s/output" % (self.base_path)
        self.report_folder = "%s/reports_and_stats" % (self.output_folder)
        self.raw_stat_data_folder = "%s/stats_data_json" % (self.output_folder)
        self._set_input_folder_names()
        self._set_read_mapping_folder_names()
        self._set_coverage_folder_names()
        self._set_gene_quanti_folder_names()

    def _set_input_folder_names(self):
        self.read_fasta_folder = "%s/reads" % self.input_folder
        self.genome_folder = "%s/genomes" % self.input_folder
        self.annotation_folder = "%s/annotation_files" % self.input_folder

    def _set_read_mapping_folder_names(self):
        self.read_mapping_index_folder = "%s/read_mappings-index" % (
            self.output_folder)
        self.read_mappings_folder = "%s/read_mappings-mappings" % (
            self.output_folder)
        self.processed_reads_folder = "%s/read_mappings-processed_reads" % (
            self.output_folder)
        self.unmapped_reads_folder = "%s/read_mappings-unmapped_reads" % (
            self.output_folder)

    def _set_coverage_folder_names(self):
        self.coverage_folder = "%s/coverages" % self.output_folder
        self.coverage_folder_norm_reads = "%s/coverages-read_normalized" % (
            self.output_folder)
        self.coverage_folder_norm_reads_mil = (
            "%s/coverages-read_normalized_times_million" % self.output_folder)

    def _set_gene_quanti_folder_names(self):
        self.gene_quanti_folder = (
            "%s/gene_wise_quantifications" % self.output_folder)
        self.gene_wise_quanti_combined_path = (
            "%s/Gene_wise_quantifications_combined.csv" % self.output_folder)

    def _set_static_file_names(self):
        """Set name of common files."""
        self.config_file = "%s/rapl.config" % self.base_path
        self.read_processing_stats = "%s/read_processing.json" % (
            self.raw_stat_data_folder)
        self.read_file_stats = "%s/read_file_stats.txt" % (
            self.report_folder)
        self.genome_file_stats = "%s/genome_file_stats.txt" % (
            self.report_folder)
        self.annotation_file_stats = "%s/annotation_file_stats.txt" % (
            self.report_folder)
        self.read_mapping_stat_file = "%s/read_mapping_stats.csv" % (
            self.report_folder)
        self.report_tex_file = "%s/report.tex" % (
            self.report_folder)
        self.lib_genome_read_mapping_summary = (
            "%s/mapping_and_mapped_read_coutings_per_genome_file.csv" % (
                self.report_folder))
        self.mapping_length_hist_r_file = (
            "%s/mapping_length_distributions.R" % self.report_folder)
        self.mapping_length_hist_pdf_file = (
            "%s/mapping_length_distributions.pdf" % self.report_folder)
        self.index_file_path = "%s/index.idx" % self.read_mapping_index_folder

    def _get_sorted_folder_content(self, folder):
        """Return the sorted file list of a folder"""
        return(list(filter(lambda file_name: not file_name.endswith("~"),
                      sorted(os.listdir(folder)))))

    def _get_read_file_names(self):
        """Read the names of the read files."""
        return(self._get_sorted_folder_content(self.read_fasta_folder))

    def _get_genome_file_names(self):
        """Read the names of genome files."""
        return(self._get_sorted_folder_content(self.genome_folder))

    def _get_annotation_file_names(self):
        """Read the names of annotation files."""
        return(self._get_sorted_folder_content(self.annotation_folder))

    def required_folders(self):
        return(self._required_base_folders() +
               self._required_input_folders() +
               self._required_read_mapping_folders() +
               self._required_coverage_folders() +
               self._required_gene_quanti_folder())

    def _required_base_folders(self):
        return([self.input_folder, self.output_folder, self.report_folder,
                self.raw_stat_data_folder])

    def _required_input_folders(self):
        return([self.read_fasta_folder, self.genome_folder,
                self.annotation_folder])

    def _required_read_mapping_folders(self):
        return([self.read_mappings_folder, self.processed_reads_folder,
                self.unmapped_reads_folder, self.read_mapping_index_folder])

    def _required_coverage_folders(self):
        return([self.coverage_folder, self.coverage_folder_norm_reads,
                self.coverage_folder_norm_reads_mil])

    def _required_gene_quanti_folder(self):
        return([self.gene_quanti_folder])

    def set_read_files_dep_file_lists(self, read_files):
        self.read_file_paths = self._path_list(self.read_fasta_folder, read_files)
        self.processed_read_file_paths = self._path_list(
            self.processed_reads_folder, read_files, appendix="_processed.fa")
        self.read_mapping_result_sam_paths = self._path_list(
            self.read_mappings_folder, read_files, appendix="_mappings.sam")
        # samtool appends ".bam" so only the prefix is required
        self.read_mapping_result_bam_prefixes_paths = self._path_list(
            self.read_mappings_folder, read_files, appendix="_mappings")
        self.read_mapping_result_bam_paths = self._path_list(
            self.read_mappings_folder, read_files, appendix="_mappings.bam")
        self.unmapped_reads_paths = self._path_list(
            self.unmapped_reads_folder, read_files, appendix="_unmapped.fa")

    def set_genome_paths(self, genome_files):
        self.genome_file_paths = self._path_list(
            self.genome_folder, genome_files)

    def set_annotation_paths(self, annotation_files):
        self.annotation_file_paths = self._path_list(
            self.annotation_folder, annotation_files)

    def _path_list(self, folder, files, appendix=""):
        return(["%s/%s%s" % (folder, file, appendix) for file in files])

    def gene_quanti_path(self, read_file_name, annotation_file):
        return("%s/%s_to_%s.csv" % (
            self.gene_quanti_folder, read_file_name, annotation_file))
