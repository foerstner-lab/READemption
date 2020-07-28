import os
import sys


class Paths(object):
    def __init__(self, base_path="."):
        self.base_path = base_path
        self._set_folder_names()
        self._set_static_files()

    def _set_folder_names(self):
        """Set the name of folders used in a project."""
        self.input_folder = f"{self.base_path}/input"
        self.output_folder = f"{self.base_path}/output"
        self._set_input_folder_names()
        self._set_read_alignment_folder_names()
        self._set_coverage_folder_names()
        self._set_gene_quanti_folder_names()
        self._set_deseq_folder_names()
        self._set_viz_align_folder_names()
        self._set_viz_gene_quanti_folder_names()
        self._set_viz_deseq_folder_names()

    def _set_input_folder_names(self):
        self.read_fasta_folder = f"{self.input_folder}/reads"
        self.ref_seq_folder = f"{self.input_folder}/reference_sequences"
        self.annotation_folder = f"{self.input_folder}/annotations"

    def _set_read_alignment_folder_names(self):
        self.align_base_folder = f"{self.output_folder}/align"
        self.read_alignment_index_folder = f"{self.align_base_folder}/index"
        self.read_alignments_folder = f"{self.align_base_folder}/alignments"
        self.processed_reads_folder = (
            f"{self.align_base_folder}/processed_reads"
        )
        self.unaligned_reads_folder = (
            f"{self.align_base_folder}/unaligned_reads"
        )
        self.align_report_folder = f"{self.align_base_folder}/reports_and_stats"
        self.raw_stat_data_folder = (
            f"{self.align_report_folder}/stats_data_json"
        )

    def _set_coverage_folder_names(self):
        self.coverage_base_folder = f"{self.output_folder}/coverage"
        self.coverage_raw_folder = f"{self.coverage_base_folder}/coverage-raw"
        self.coverage_tnoar_min_norm_folder = (
            f"{self.coverage_base_folder}/coverage-tnoar_min_normalized"
        )
        self.coverage_tnoar_mil_norm_folder = (
            f"{self.coverage_base_folder}/coverage-tnoar_mil_normalized"
        )

    def _set_gene_quanti_folder_names(self):
        self.gene_quanti_base_folder = f"{self.output_folder}/gene_quanti"
        self.gene_quanti_per_lib_folder = (
            f"{self.gene_quanti_base_folder}/gene_quanti_per_lib"
        )
        self.gene_quanti_combined_folder = (
            f"{self.gene_quanti_base_folder}/gene_quanti_combined"
        )
        self.gene_wise_quanti_combined_path = f"{self.gene_quanti_combined_folder}/gene_wise_quantifications_combined.csv"
        self.gene_wise_quanti_combined_rpkm_path = f"{self.gene_quanti_combined_folder}/gene_wise_quantifications_combined_rpkm.csv"
        self.gene_wise_quanti_combined_tnoar_path = f"{self.gene_quanti_combined_folder}/gene_wise_quantifications_combined_tnoar.csv"
        self.gene_wise_quanti_combined_tpm_path = f"{self.gene_quanti_combined_folder}/gene_wise_quantifications_combined_tpm.csv"

    def _set_deseq_folder_names(self):
        self.deseq_base_folder = f"{self.output_folder}/deseq"
        self.deseq_raw_folder = f"{self.deseq_base_folder}/deseq_raw"
        self.deseq_extended_folder = (
            f"{self.deseq_base_folder}/deseq_with_annotations"
        )

    def _set_viz_align_folder_names(self):
        self.viz_align_base_folder = f"{self.output_folder}/viz_align"
        self.viz_align_input_read_length_plot_path = (
            f"{self.viz_align_base_folder}/input_reads_length_distributions.pdf"
        )
        self.viz_align_processed_reads_length_plot_path = f"{self.viz_align_base_folder}/processed_reads_length_distributions.pdf"

    def _set_viz_gene_quanti_folder_names(self):
        self.viz_gene_quanti_base_folder = (
            f"{self.output_folder}/viz_gene_quanti"
        )
        self.viz_gene_quanti_scatter_plot_path = (
            f"{self.viz_gene_quanti_base_folder}/expression_scatter_plots.pdf"
        )
        self.viz_gene_quanti_rna_classes_plot_path = (
            f"{self.viz_gene_quanti_base_folder}/rna_class_sizes.pdf"
        )

    def _set_viz_deseq_folder_names(self):
        self.viz_deseq_base_folder = f"{self.output_folder}/viz_deseq"
        self.viz_deseq_scatter_plot_path = (
            f"{self.viz_deseq_base_folder}/MA_plots.pdf"
        )
        self.viz_deseq_volcano_plot_path = f"{self.viz_deseq_base_folder}/volcano_plots_log2_fold_change_vs_p-value.pdf"
        self.viz_deseq_volcano_plot_adj_path = f"{self.viz_deseq_base_folder}/volcano_plots_log2_fold_change_vs_adjusted_p-value.pdf"

    def _set_static_files(self):
        """Set name of common files."""
        self.read_processing_stats_path = (
            f"{self.raw_stat_data_folder}/read_processing.json"
        )
        self.primary_read_aligner_stats_path = (
            f"{self.raw_stat_data_folder}/read_alignments_primary_aligner.json"
        )

        self.read_realigner_stats_path = (
            f"{self.raw_stat_data_folder}/read_alignments_realigner.json"
        )
        self.read_alignments_stats_path = (
            f"{self.raw_stat_data_folder}/read_alignments_final.json"
        )
        self.read_file_stats = (
            f"{self.align_report_folder}/input_read_stats.txt"
        )
        self.ref_seq_file_stats = (
            f"{self.align_report_folder}/reference_sequences_stats.txt"
        )
        self.read_alignment_stats_table_path = (
            f"{self.align_report_folder}/read_alignment_stats.csv"
        )
        self.index_path = f"{self.read_alignment_index_folder}/index.idx"
        self.deseq_script_path = f"{self.deseq_raw_folder}/deseq.R"
        self.deseq_pca_heatmap_path = (
            f"{self.deseq_raw_folder}/sample_comparison_pca_heatmap.pdf"
        )
        self.deseq_tmp_session_info_script = f"{self.deseq_raw_folder}/tmp.R"
        self.deseq_session_info = f"{self.deseq_raw_folder}/R_session_info.txt"
        self.version_path = f"{self.align_report_folder}/version_log.txt"

    def _get_sorted_folder_content(self, folder):
        """Return the sorted file list of a folder"""
        return list(
            filter(
                lambda file: not (
                    file.endswith("~") or os.path.basename(file).startswith(".")
                ),
                sorted(os.listdir(folder)),
            )
        )

    def get_read_files(self):
        """Read the names of the read files."""
        return self._get_sorted_folder_content(self.read_fasta_folder)

    def get_read_file_pairs(self):
        """Read the names of the read files as paired for paired end reads."""
        read_files = self._get_sorted_folder_content(self.read_fasta_folder)
        if len(read_files) % 2 != 0:
            sys.stderr.write(
                "Error: Number of files is unequal. This cannot be "
                "the case for paired end run data.\n"
            )
            sys.exit(1)
        for read_file in read_files:
            if not ("_p1" in read_file or "_p2" in read_file):
                sys.stderr.write(
                    "Error: You specified this as paired end run data but the "
                    "file name '{}' does not contain '_p1' or '_p2'.\n"
                )
                sys.exit(1)
        read_file_pairs = list(
            read_file_pair
            for read_file_pair in zip(read_files[::2], read_files[1::2])
        )
        for read_file_pair in read_file_pairs:
            if read_file_pair[0] != read_file_pair[1].replace("_p2.", "_p1."):
                sys.stderr.write(
                    "Error: Cannot detect pairs of paired end "
                    "sequenceing files. Please check that file "
                    "names contain '_p1' and '_p2'.\n"
                )
                sys.exit(1)
        return read_file_pairs

    def get_lib_names_single_end(self):
        """Extract the suffux free name of single end libraries"""
        return [
            self._clean_file_name(file_name)
            for file_name in self.get_read_files()
        ]

    def get_lib_names_paired_end(self):
        """Extract the suffux free name of paired end libraries

        For each pair of read files one library name is returned.

        """
        p1_names = [
            self._clean_file_name(file_name)
            for file_name in self.get_read_files()
        ][::2]
        for p1_name in p1_names:
            if not p1_name.endswith("_p1"):
                sys.stderr.write(
                    f"Error: File '{p1_name}' should end with '_p1' but "
                    "does not. Please check file name convention "
                    "for paired end reads.\n"
                )
                sys.exit(1)
        return [p1_name[:-3] for p1_name in p1_names]

    def _clean_file_name(self, file_name):
        for suffix in [
            "bz2",
            "BZ",
            "gz",
            "GZ",
            "fa",
            "fasta",
            "FA",
            "FASTA",
            "fastq",
            "fq",
            "FQ",
            "FASTQ",
        ]:
            if file_name.endswith(suffix):
                suffix = "." + suffix
                file_name = file_name[: -len(suffix)]
        return file_name

    def get_ref_seq_files(self):
        """Read the names of reference sequence files."""
        return self._get_sorted_folder_content(self.ref_seq_folder)

    def get_annotation_files(self):
        """Read the names of annotation files."""
        return self._get_sorted_folder_content(self.annotation_folder)

    def required_folders(self):
        return (
            self.required_base_folders()
            + self.required_input_folders()
            + self.required_read_alignment_folders()
            + self.required_coverage_folders()
            + self.required_gene_quanti_folders()
            + self.required_deseq_folders()
            + self.required_viz_align_folders()
            + self.required_viz_gene_quanti_folders()
            + self.required_viz_deseq_folders()
        )

    def required_base_folders(self):
        return [self.input_folder, self.output_folder]

    def required_input_folders(self):
        return [
            self.read_fasta_folder,
            self.ref_seq_folder,
            self.annotation_folder,
        ]

    def required_read_alignment_folders(self):
        return [
            self.align_base_folder,
            self.read_alignments_folder,
            self.processed_reads_folder,
            self.unaligned_reads_folder,
            self.read_alignment_index_folder,
            self.align_report_folder,
            self.raw_stat_data_folder,
        ]

    def required_coverage_folders(self):
        return [
            self.coverage_base_folder,
            self.coverage_raw_folder,
            self.coverage_tnoar_min_norm_folder,
            self.coverage_tnoar_mil_norm_folder,
        ]

    def required_gene_quanti_folders(self):
        return [
            self.gene_quanti_base_folder,
            self.gene_quanti_per_lib_folder,
            self.gene_quanti_combined_folder,
        ]

    def required_deseq_folders(self):
        return [
            self.deseq_base_folder,
            self.deseq_raw_folder,
            self.deseq_extended_folder,
        ]

    def required_viz_align_folders(self):
        return [self.viz_align_base_folder]

    def required_viz_gene_quanti_folders(self):
        return [self.viz_gene_quanti_base_folder]

    def required_viz_deseq_folders(self):
        return [self.viz_deseq_base_folder]

    def set_read_files_dep_file_lists_single_end(self, read_files, lib_names):
        self.read_paths = self._path_list(self.read_fasta_folder, read_files)
        self.processed_read_paths = self._path_list(
            self.processed_reads_folder, lib_names, appendix="_processed.fa.gz"
        )
        self.unaligned_reads_paths = self._path_list(
            self.unaligned_reads_folder, lib_names, appendix="_unaligned.fa"
        )
        self.realigned_unaligned_reads_paths = self._path_list(
            self.unaligned_reads_folder,
            lib_names,
            appendix="_unaligned_after_realignment.fa",
        )
        self._set_alignment_paths(lib_names)

    def set_read_files_dep_file_lists_paired_end(
        self, read_file_pairs, lib_names
    ):
        self.read_path_pairs = [
            self._path_list(self.read_fasta_folder, read_file_pair)
            for read_file_pair in read_file_pairs
        ]
        self.processed_read_path_pairs = [
            self._path_list(
                self.processed_reads_folder,
                [
                    self._clean_file_name(read_file)
                    for read_file in read_file_pair
                ],
                appendix="_processed.fa.gz",
            )
            for read_file_pair in read_file_pairs
        ]
        # The read of both files that are not matchend will be dumped
        # together into on file. Due to this there is only one file
        # per pair.
        self.unaligned_reads_paths = self._path_list(
            self.unaligned_reads_folder, lib_names, appendix="_unaligned.fa"
        )
        self.realigned_unaligned_reads_paths = self._path_list(
            self.unaligned_reads_folder,
            lib_names,
            appendix="_unaligned_after_realignment.fa",
        )
        self._set_alignment_paths(lib_names)

    def _set_alignment_paths(self, lib_names):
        ###
        # For the primary mapper
        # self.primary_read_aligner_sam_paths = self._path_list(
        #    self.read_alignments_folder, lib_names,
        #    appendix="_alignments_primary_aligner.sam")
        self.primary_read_aligner_bam_paths = self._path_list(
            self.read_alignments_folder,
            lib_names,
            appendix="_alignments_primary_aligner.bam",
        )
        # samtool appends ".bam" so only the prefix is required
        # self.primary_read_aligner_bam_prefix_paths = self._path_list(
        #    self.read_alignments_folder,
        #    lib_names, appendix="_alignments_primary_aligner")
        ###
        # For the remapper
        # self.read_realigner_tmp_sam_paths = self._path_list(
        #    self.read_alignments_folder, lib_names,
        #    appendix="_alignments_tmp_for_realigner.sam")
        # self.read_realigner_sam_paths = self._path_list(
        #    self.read_alignments_folder, lib_names,
        #    appendix="_alignments_realigner.sam")
        self.read_realigner_bam_paths = self._path_list(
            self.read_alignments_folder,
            lib_names,
            appendix="_alignments_realigner.bam",
        )
        # self.read_realigner_bam_prefixes_paths = self._path_list(
        #    self.read_alignments_folder,
        #    lib_names, appendix="_alignments_realigner")
        ###
        # For the cross-aligned cleaned
        self.read_alignment_bam_cross_cleaned_tmp_paths = self._path_list(
            self.read_alignments_folder,
            lib_names,
            appendix="_alignments_tmp_crossmapped_cleaned.bam",
        )
        self.read_alignment_bam_with_crossmappings_paths = self._path_list(
            self.read_alignments_folder,
            lib_names,
            appendix="_alignments_potententially_with_crossmappings.bam",
        )
        self.crossmapped_reads_paths = self._path_list(
            self.read_alignments_folder,
            lib_names,
            appendix="_crossmapped_reads.txt",
        )
        ###
        # For the final (merged) version
        self.read_alignment_bam_paths = self._path_list(
            self.read_alignments_folder,
            lib_names,
            appendix="_alignments_final.bam",
        )
        self.read_alignment_bam_prefix_paths = self._path_list(
            self.read_alignments_folder, lib_names, appendix="_alignments_final"
        )

    def set_ref_seq_paths(self, ref_seq_files):
        self.ref_seq_paths = self._path_list(self.ref_seq_folder, ref_seq_files)

    def set_annotation_paths(self, annotation_files):
        self.annotation_paths = self._path_list(
            self.annotation_folder, annotation_files
        )

    def _path_list(self, folder, files, appendix=""):
        return [f"{folder}/{file}{appendix}" for file in files]

    def gene_quanti_path(self, read_file, annotation_file):
        return f"{self.gene_quanti_per_lib_folder}/{read_file}_to_{annotation_file}.csv"

    def wiggle_file_raw_path(self, read_file, strand, multi=None, div=None):
        return self._wiggle_file_path(
            self.coverage_raw_folder, read_file, strand, multi=None, div=None
        )

    def wiggle_file_tnoar_norm_min_path(
        self, read_file, strand, multi=None, div=None
    ):
        return self._wiggle_file_path(
            self.coverage_tnoar_min_norm_folder, read_file, strand, multi, div
        )

    def wiggle_file_tnoar_norm_mil_path(
        self, read_file, strand, multi=None, div=None
    ):
        return self._wiggle_file_path(
            self.coverage_tnoar_mil_norm_folder, read_file, strand, multi, div
        )

    def _wiggle_file_path(
        self, folder, read_file, strand, multi=None, div=None
    ):
        path = f"{folder}/{read_file}"
        if not div is None:
            path += f"_div_by_{div}.1f"
        if not multi is None:
            path += f"_multi_by_{multi}.1f"
        path += f"_{strand}.wig"
        return path
