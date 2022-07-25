from collections import defaultdict
import json
import os
import sys


class PathCreator:
    def __init__(self, base_path, species_information):
        self.ref_seq_paths_by_species = None
        self.base_path = base_path
        self.config_file = f"{self.base_path}/config.json"
        self.species_information = species_information
        # The species information is either taken from a command line argument
        # or if that argument was not given ("None"), the speces information
        # is retrieved from the config file
        if self.species_information:
            self.species_folder_prefixes_and_display_names = (
                self._get_species_from_cml_argument(self.species_information)
            )
        else:
            self.species_folder_prefixes_and_display_names = (
                self._get_species_from_config(self.config_file)
            )
        (
            self.species_folder_prefixes,
            self.species_display_names,
        ) = self._set_species_folder_prefixes_and_display_names(
            self.species_folder_prefixes_and_display_names
        )
        self.prefix_folder_name_connector = "_"
        self._set_folder_names()
        self._set_static_files()

    def _get_species_from_cml_argument(self, species_information: str) -> dict:
        """
        :param species_information: The argument regarding the species
        information that was provided by the user
        :return: Dictionary containing the species information
        """
        species_folder_prefixes_and_display_names = {}
        for sp in species_information:
            species_folder_prefix, display_name = sp.split("=")
            species_folder_prefixes_and_display_names[
                species_folder_prefix
            ] = display_name
        return species_folder_prefixes_and_display_names

    def _get_species_from_config(self, config_file: str) -> dict:
        """
        :param config_file: The path of the config file
        :return: Dictionary containing the species information
        """
        with open(config_file, "r") as config_file:
            config = json.loads(config_file.read())
            return config["species"]

    def _set_species_folder_prefixes_and_display_names(
        self, species_folder_prefixes_and_display_names
    ):
        species_folder_prefixes = (
            species_folder_prefixes_and_display_names.keys()
        )
        species_display_names = (
            species_folder_prefixes_and_display_names.values()
        )
        return species_folder_prefixes, species_display_names

    def _set_folder_names(self):
        """Set the name of folders used in a project."""
        self.input_folder = f"{self.base_path}/input"
        self.output_folder = f"{self.base_path}/output"
        self._set_input_folder_names()
        self._set_read_alignment_folder_names()
        # self._set_coverage_folder_names() #TODO needs to be set in controller
        self._set_gene_quanti_folder_names()
        self._set_deseq_folders_and_files_by_species()
        self._set_viz_align_folder_names()
        self._set_viz_gene_quanti_folder_names()
        self._set_viz_deseq_folder_names()

    def _set_input_folder_names(self):
        self.input_folder
        self.read_fasta_folder = f"{self.input_folder}/reads"
        self._set_ref_seq_folders()
        self._set_annotation_folders()

    def _set_ref_seq_folders(self):
        self.ref_seq_folders_by_species = {}
        for prefix in self.species_folder_prefixes:
            prefix_and_connector = prefix + self.prefix_folder_name_connector
            if len(self.species_folder_prefixes) and prefix == " ":
                prefix_and_connector = ""
            self.ref_seq_folders_by_species[
                prefix
            ] = f"{self.input_folder}/{prefix_and_connector}reference_sequences"

    def _set_annotation_folders(self):
        self.annotation_folders_by_species = {}
        for prefix in self.species_folder_prefixes:
            prefix_and_connector = prefix + self.prefix_folder_name_connector
            if len(self.species_folder_prefixes) and prefix == " ":
                prefix_and_connector = ""
            self.annotation_folders_by_species[
                prefix
            ] = f"{self.input_folder}/{prefix_and_connector}annotations"

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

    def set_coverage_folder_and_file_names(
        self,
        strands,
        lib_names,
        read_files_aligned_read_freq_and_min_reads_aligned_by_species,
    ):
        self.coverage_folders_by_species = {}
        self.coverage_files_by_species = {}
        for prefix in self.species_folder_prefixes:
            prefix_and_connector = prefix + self.prefix_folder_name_connector
            if len(self.species_folder_prefixes) and prefix == " ":
                prefix_and_connector = ""
            coverage_species_folders = {}
            coverage_species_folders[
                "coverage_raw_folder"
            ] = f"{self.output_folder}/{prefix_and_connector}coverage-raw"
            coverage_species_folders[
                "coverage_tnoar_min_norm_folder"
            ] = f"{self.output_folder}/{prefix_and_connector}coverage-tnoar_min_normalized"
            coverage_species_folders[
                "coverage_tnoar_mil_norm_folder"
            ] = f"{self.output_folder}/{prefix_and_connector}coverage-tnoar_mil_normalized"
            self.coverage_folders_by_species[prefix] = coverage_species_folders

            coverage_species_files = {}
            for lib in lib_names:
                # Add wiggle raw paths
                coverage_species_files[lib] = {}
                coverage_species_files[lib]["wiggle_file_raw_path"] = {}
                for strand in strands:
                    coverage_species_files[lib]["wiggle_file_raw_path"][
                        strand
                    ] = self._wiggle_file_path(
                        coverage_species_folders["coverage_raw_folder"],
                        lib,
                        strand,
                    )

                # Add wiggle tnoar norm min path tnoar_norm_min
                coverage_species_files[lib][
                    "wiggle_file_tnoar_norm_min_path"
                ] = {}
                for strand in strands:
                    coverage_species_files[lib][
                        "wiggle_file_tnoar_norm_min_path"
                    ][strand] = self._wiggle_file_path(
                        coverage_species_folders[
                            "coverage_tnoar_min_norm_folder"
                        ],
                        lib,
                        strand,
                        multi=read_files_aligned_read_freq_and_min_reads_aligned_by_species[
                            prefix
                        ][
                            "min_no_of_aligned_reads"
                        ],
                        div=read_files_aligned_read_freq_and_min_reads_aligned_by_species[
                            prefix
                        ][
                            "read_files_aligned_read_freq"
                        ][
                            lib
                        ],
                    )
                # Add wiggle tnoar norm min path tnoar_norm_mil
                coverage_species_files[lib][
                    "wiggle_file_tnoar_norm_mil_path"
                ] = {}
                for strand in strands:
                    coverage_species_files[lib][
                        "wiggle_file_tnoar_norm_mil_path"
                    ][strand] = self._wiggle_file_path(
                        coverage_species_folders[
                            "coverage_tnoar_mil_norm_folder"
                        ],
                        lib,
                        strand,
                        multi=1000000,
                        div=read_files_aligned_read_freq_and_min_reads_aligned_by_species[
                            prefix
                        ][
                            "read_files_aligned_read_freq"
                        ][
                            lib
                        ],
                    )

            self.coverage_files_by_species[prefix] = coverage_species_files

            # weiter
            # def wiggle_file_raw_path(self, read_file, strand, multi=None, div=None):
            #    return self._wiggle_file_path(
            #        self.coverage_raw_folder, read_file, strand, multi=None, div=None
            #    )

            #    def wiggle_file_raw_path(self, read_file, strand, multi=None, div=None):

    #        return self._wiggle_file_path(
    #            self.coverage_raw_folder, read_file, strand, multi=None, div=None
    #        )

    def _wiggle_file_path(
        self, folder, read_file, strand, multi=None, div=None
    ):
        path = f"{folder}/{read_file}"
        if not div is None:
            path += f"_div_by_{div:.1f}"
        if not multi is None:
            path += f"_multi_by_{multi:.1f}"
        path += f"_{strand}.wig"
        return path

    def _set_gene_quanti_folder_names(self):
        self.gene_quanti_folders_by_species = {}
        self.gene_quanti_files_by_species = {}
        for prefix in self.species_folder_prefixes:
            prefix_and_connector = prefix + self.prefix_folder_name_connector
            if len(self.species_folder_prefixes) and prefix == " ":
                prefix_and_connector = ""
            gene_quanti_species_folders = {}
            gene_quanti_species_folders[
                "gene_quanti_base_folder"
            ] = f"{self.output_folder}/{prefix_and_connector}gene_quanti"
            gene_quanti_species_folders[
                "gene_quanti_per_lib_folder"
            ] = f"{self.output_folder}/{prefix_and_connector}gene_quanti_per_lib"
            gene_quanti_species_folders[
                "gene_quanti_combined_folder"
            ] = f"{self.output_folder}/{prefix_and_connector}gene_quanti_combined"
            self.gene_quanti_folders_by_species[
                prefix
            ] = gene_quanti_species_folders

            gene_quanti_species_files = {}
            gene_quanti_species_files[
                "gene_wise_quanti_combined_path"
            ] = f"{gene_quanti_species_folders['gene_quanti_combined_folder']}/gene_wise_quantifications_combined.csv"
            gene_quanti_species_files[
                "gene_wise_quanti_combined_rpkm_path"
            ] = f"{gene_quanti_species_folders['gene_quanti_combined_folder']}/gene_wise_quantifications_combined_rpkm.csv"
            gene_quanti_species_files[
                "gene_wise_quanti_combined_tnoar_path"
            ] = f"{gene_quanti_species_folders['gene_quanti_combined_folder']}/gene_wise_quantifications_combined_tnoar.csv"
            gene_quanti_species_files[
                "gene_wise_quanti_combined_tpm_path"
            ] = f"{gene_quanti_species_folders['gene_quanti_combined_folder']}/gene_wise_quantifications_combined_tpm.csv"
            self.gene_quanti_files_by_species[
                prefix
            ] = gene_quanti_species_files

    def _set_deseq_folders_and_files_by_species(self):
        self.deseq_folders_by_species = {}
        self.deseq_files_by_species = {}
        for prefix in self.species_folder_prefixes:
            prefix_and_connector = prefix + self.prefix_folder_name_connector
            if len(self.species_folder_prefixes) and prefix == " ":
                prefix_and_connector = ""
            deseq_species_folders = {}
            deseq_species_folders[
                "deseq_base_folder"
            ] = f"{self.output_folder}/{prefix_and_connector}deseq"
            deseq_species_folders[
                "deseq_raw_folder"
            ] = f"{self.output_folder}/{prefix_and_connector}deseq/deseq_raw"
            deseq_species_folders[
                "deseq_extended_folder"
            ] = f"{self.output_folder}/{prefix_and_connector}deseq/deseq_with_annotations"

            self.deseq_folders_by_species[prefix] = deseq_species_folders

            deseq_species_files = {}
            deseq_species_files[
                "deseq_script_path"
            ] = f"{deseq_species_folders['deseq_raw_folder']}/deseq.R"

            deseq_species_files[
                "deseq_pca_heatmap_path"
            ] = f"{deseq_species_folders['deseq_raw_folder']}/sample_comparison_pca_heatmap.pdf"

            deseq_species_files[
                "deseq_tmp_session_info_script"
            ] = f"{deseq_species_folders['deseq_raw_folder']}/tmp.R"

            deseq_species_files[
                "deseq_session_info"
            ] = f"{deseq_species_folders['deseq_raw_folder']}/R_session_info.txt"

            deseq_species_files[
                "deseq_session_info"
            ] = f"{deseq_species_folders['deseq_raw_folder']}/R_session_info.txt"

            self.deseq_files_by_species[prefix] = deseq_species_files

    def _set_viz_align_folder_names(self):
        self.viz_align_read_lengths_folder = (
            f"{self.output_folder}/read_lengths_viz_align"
        )
        self.viz_align_all_folder = (
            f"{self.output_folder}/all_species_viz_align"
        )
        self.viz_align_input_read_length_plot_path = f"{self.viz_align_read_lengths_folder}/input_reads_length_distributions.pdf"
        self.viz_align_processed_reads_length_plot_path = f"{self.viz_align_read_lengths_folder}/processed_reads_length_distributions.pdf"
        self._set_viz_align_folders_by_species()

    def _set_viz_align_folders_by_species(self):
        self.viz_align_folders_by_species = {}
        self.viz_align_aligned_reads_by_species_paths = {}
        for prefix in self.species_folder_prefixes:
            prefix_and_connector = prefix + self.prefix_folder_name_connector
            if len(self.species_folder_prefixes) and prefix == " ":
                prefix_and_connector = ""
            self.viz_align_folders_by_species[
                prefix
            ] = f"{self.output_folder}/{prefix_and_connector}viz_align"

            self.viz_align_aligned_reads_by_species_paths[
                prefix
            ] = f"{self.output_folder}/{prefix_and_connector}viz_align/{prefix_and_connector}aligned_reads"

    def get_viz_align_folders_by_species(self):
        self._set_viz_align_folders_by_species()
        return self.viz_align_folders_by_species

    def _set_viz_gene_quanti_folder_names(self):
        self.viz_gene_quanti_folders_by_species = {}
        self.viz_gene_quanti_files_by_species = {}
        for prefix in self.species_folder_prefixes:
            prefix_and_connector = prefix + self.prefix_folder_name_connector
            if len(self.species_folder_prefixes) and prefix == " ":
                prefix_and_connector = ""
            viz_gene_quanti_species_folders = {}
            viz_gene_quanti_species_folders[
                "viz_gene_quanti_base_folder"
            ] = f"{self.output_folder}/{prefix_and_connector}viz_gene_quanti"
            self.viz_gene_quanti_folders_by_species[
                prefix
            ] = viz_gene_quanti_species_folders

            viz_gene_quanti_species_files = {}
            viz_gene_quanti_species_files[
                "viz_gene_quanti_scatter_plot_path"
            ] = f"{viz_gene_quanti_species_folders['viz_gene_quanti_base_folder']}/expression_scatter_plots.pdf"
            viz_gene_quanti_species_files[
                "rna_classes_plot_path"
            ] = f"{viz_gene_quanti_species_folders['viz_gene_quanti_base_folder']}/rna_class_sizes.pdf"

            self.viz_gene_quanti_files_by_species[
                prefix
            ] = viz_gene_quanti_species_files

    def _set_viz_deseq_folder_names(self):
        self.viz_deseq_folders_by_species = {}
        self.viz_deseq_files_by_species = {}
        for prefix in self.species_folder_prefixes:
            prefix_and_connector = prefix + self.prefix_folder_name_connector
            if len(self.species_folder_prefixes) and prefix == " ":
                prefix_and_connector = ""
            viz_deseq_species_folders = {}
            viz_deseq_species_folders[
                "viz_deseq_base_folder"
            ] = f"{self.output_folder}/{prefix_and_connector}viz_deseq"
            self.viz_deseq_folders_by_species[
                prefix
            ] = viz_deseq_species_folders

            viz_deseq_species_files = {}
            viz_deseq_species_files[
                "viz_deseq_scatter_plot_path"
            ] = f"{viz_deseq_species_folders['viz_deseq_base_folder']}/MA_plots.pdf"
            viz_deseq_species_files[
                "viz_deseq_volcano_plot_path"
            ] = f"{viz_deseq_species_folders['viz_deseq_base_folder']}/volcano_plots_log2_fold_change_vs_p-value.pdf"
            viz_deseq_species_files[
                "viz_deseq_volcano_plot_adj_path"
            ] = f"{viz_deseq_species_folders['viz_deseq_base_folder']}/volcano_plots_log2_fold_change_vs_adjusted_p-value.pdf"

            self.viz_deseq_files_by_species[prefix] = viz_deseq_species_files

        # self.viz_deseq_base_folder = f"{self.output_folder}/viz_deseq"
        # self.viz_deseq_scatter_plot_path = (
        #    f"{self.viz_deseq_base_folder}/MA_plots.pdf"
        # )
        # self.viz_deseq_volcano_plot_path = f"{self.viz_deseq_base_folder}/volcano_plots_log2_fold_change_vs_p-value.pdf"
        # self.viz_deseq_volcano_plot_adj_path = f"{self.viz_deseq_base_folder}/volcano_plots_log2_fold_change_vs_adjusted_p-value.pdf"

    def _set_static_files(self):
        """Set name of common files."""
        self.read_processing_stats_path = (
            f"{self.raw_stat_data_folder}/read_processing.json"
        )
        self.read_alignments_stats_path = (
            f"{self.raw_stat_data_folder}/read_alignments_final.json"
        )
        self.fragment_alignments_stats_path = (
            f"{self.raw_stat_data_folder}/fragment_alignments_final.json"
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
        self.read_alignment_stats_table_transposed_path = (
            f"{self.align_report_folder}/read_alignment_stats_transposed.csv"
        )
        self.fragment_alignment_stats_table_path = (
            f"{self.align_report_folder}/fragment_alignment_stats.csv"
        )
        self.fragment_alignment_stats_table_transposed_path = (
            f"{self.align_report_folder}/fragment_alignment_stats_transposed.csv"
        )
        self.index_path = f"{self.read_alignment_index_folder}/index.idx"
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
        """Extract the suffix free name of single end libraries"""
        return [
            self._clean_file_name(file_name)
            for file_name in self.get_read_files()
        ]

    def get_lib_names_paired_end(self):
        """Extract the suffix free name of paired end libraries

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

    def set_ref_seq_paths_by_species(self) -> None:
        """
        sets the attribute ref_seq_paths_by_species that is a dictionary
        containing the reference sequence paths sorted by species.
        E.g.:
        {'human':
        ['reademption_analysis_dual/input/human_reference_sequences/human.fa',
        'reademption_analysis_dual/input/human_reference_sequences/human_contig.fa],

        'e_coli':
        ['reademption_analysis_dual/input/e_coli_reference_sequences/e_coli.fa]}
        """
        self.ref_seq_paths_by_species = {}
        for (
            sp,
            species_ref_seq_folder,
        ) in self.ref_seq_folders_by_species.items():
            self.ref_seq_paths_by_species[sp] = self._path_list(
                species_ref_seq_folder,
                self._get_sorted_folder_content(species_ref_seq_folder),
            )

    def get_ref_seq_files(self) -> list:
        """
        extracts the names of the reference sequence files from the dictionary
        self.ref_seq_folders_by_species and writes them to a list.
        E.g.:
        ['human.fa, e_coli.fa]
        :return: a list containing all reference sequences
        """
        if not self.ref_seq_paths_by_species:
            self.set_ref_seq_paths_by_species()
        ref_seq_files = []
        for ref_seq_folder in self.ref_seq_folders_by_species.values():
            for ref_seq in self._get_sorted_folder_content(ref_seq_folder):
                ref_seq_files.append(ref_seq)
        return ref_seq_files

    def set_ref_seq_path_list(self) -> None:
        """
        sets an attribute that holds a list of all reference sequence paths.
        E.g.:
        ['reademption_analysis_dual/input/human_reference_sequences/human.fa',
        'reademption_analysis_dual/input/e_coli_reference_sequences/e_coli.fa]
        """
        if not self.ref_seq_paths_by_species:
            self.set_ref_seq_paths_by_species()
        self.ref_seq_path_list = []
        for ref_seq_paths in self.ref_seq_paths_by_species.values():
            for ref_seq_path in ref_seq_paths:
                self.ref_seq_path_list.append(ref_seq_path)

    def set_annotation_paths_by_species(self) -> None:
        """
        sets the attribute annotation_paths_by_species that is a dictionary
        containing the annotation paths sorted by species.
        E.g.:
        {'human':
        ['reademption_analysis_dual/input/human_annotations/human.gff',
        'reademption_analysis_dual/input/human_annotations/sRNA.gff],

        'e_coli':
        ['reademption_analysis_dual/input/e_coli_annotations/e_coli.gff]}
        """
        self.annotation_paths_by_species = {}
        for (
            sp,
            species_annotation_folder,
        ) in self.annotation_folders_by_species.items():
            self.annotation_paths_by_species[sp] = self._path_list(
                species_annotation_folder,
                self._get_sorted_folder_content(species_annotation_folder),
            )

    def set_annotation_files_by_species(self) -> None:
        """
        sets the attribute annotation_files_by_species that is a dictionary
        containing the annotation files sorted by species.
        E.g.:
        {'human':
            ['human.gff',
            'sRNA.gff',
        'e_coli':
            ['e_coli.gff']}
        """
        self.annotation_files_by_species = {}
        for (
            sp,
            species_annotation_folder,
        ) in self.annotation_folders_by_species.items():
            self.annotation_files_by_species[
                sp
            ] = self._get_sorted_folder_content(species_annotation_folder)

    def get_annotation_files(self) -> list:
        """
        extracts the names of the reference sequence files from the dictionary
        self.annotation_folders_by_species and writes them to a list.
        E.g.:
        ['human.fa, e_coli.fa]
        :return: a list containing all reference sequences
        """
        if not self.annotation_paths_by_species:
            self.set_annotation_paths_by_species()
        annotation_files = []
        for annotation_folder in self.annotation_folders_by_species.values():
            for annotation in self._get_sorted_folder_content(
                annotation_folder
            ):
                annotation_files.append(annotation)
        return annotation_files

    def set_annotation_path_list(self) -> None:
        """
        sets an attribute that holds a list of all reference sequence paths.
        E.g.:
        ['reademption_analysis_dual/input/human_reference_sequences/human.fa',
        'reademption_analysis_dual/input/e_coli_reference_sequences/e_coli.fa]
        """
        if not self.annotation_paths_by_species:
            self.set_annotation_paths_by_species()
        self.annotation_path_list = []
        for annotation_paths in self.annotation_paths_by_species.values():
            for annotation_path in annotation_paths:
                self.annotation_path_list.append(annotation_path)

    def required_folders(self):
        # TODO can be deleted?, because no subcommand creates all folders
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

    def required_new_project_folders(self):
        return [
            *self.required_base_folders(),
            self.read_fasta_folder,
            *self.ref_seq_folders_by_species.values(),
            *self.annotation_folders_by_species.values(),
            *self.required_read_alignment_folders(),
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
        return [*self._unpack_folder_paths(self.coverage_folders_by_species)]

    def required_viz_gene_quanti_folders(self):
        return [
            *self._unpack_folder_paths(self.viz_gene_quanti_folders_by_species)
        ]

    def required_viz_deseq_folders(self):
        return [*self._unpack_folder_paths(self.viz_deseq_folders_by_species)]

    def required_deseq_folders(self):
        return [*self._unpack_folder_paths(self.deseq_folders_by_species)]

    def _unpack_folder_paths(self, folders_by_species):
        folder_paths = []
        for species_folders in folders_by_species.values():
            for path in species_folders.values():
                folder_paths.append(path)
        return folder_paths

    def required_gene_quanti_folders(self):
        return [*self._get_gene_quanti_start_up_folders()]

    def _get_gene_quanti_start_up_folders(self):
        gene_quanti_folders = []
        for sp in self.species_folder_prefixes:
            gene_quanti_folders.append(
                self.gene_quanti_folders_by_species[sp][
                    "gene_quanti_per_lib_folder"
                ]
            )
            gene_quanti_folders.append(
                self.gene_quanti_folders_by_species[sp][
                    "gene_quanti_combined_folder"
                ]
            )
        return gene_quanti_folders

    def required_viz_align_folders(self):
        return [
            self.viz_align_read_lengths_folder,
            self.viz_align_all_folder,
            *self.viz_align_folders_by_species.values(),
        ]

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
        self.read_alignment_bam_sorted_paths = self._path_list(
            self.read_alignments_folder,
            lib_names,
            appendix="_alignments_final_sorted.bam",
        )
        self.aligned_fragments_bam_paths = self._path_list(
            self.read_alignments_folder,
            lib_names,
            appendix="_alignments_final_fragments.bam",
        )

    # TODO annotation_folder needs to be the annotation folder for the current species
    # def set_annotation_paths(self, annotation_files):
    #    self.annotation_paths = self._path_list(
    #        self.annotation_folder, annotation_files
    #    )

    def _path_list(self, folder, files, appendix=""):
        return [f"{folder}/{file}{appendix}" for file in files]

    def gene_quanti_paths_by_species(
        self, gene_quanti_per_lib_species_folder, read_file, annotation_file
    ):
        return f"{gene_quanti_per_lib_species_folder}/{read_file}_to_{annotation_file}.csv"
