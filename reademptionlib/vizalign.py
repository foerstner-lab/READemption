import json
import numpy as np
import matplotlib
import pandas as pd

import seaborn as sns
import matplotlib.patches as mpatches

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


class AlignViz(object):
    def __init__(
        self,
        lib_names,
        read_processing_stats_path,
        read_aligner_stats_path,
        read_alignment_stats_table_path,
        species_folder_prefixes_and_display_names
    ):
        self._lib_names = lib_names
        self._read_processing_stats_path = read_processing_stats_path
        self._read_aligner_stats_path = read_aligner_stats_path
        self._read_alignment_stats_table_path = read_alignment_stats_table_path
        self._species_folder_prefixes_and_display_names = species_folder_prefixes_and_display_names

    def read_stat_files(self):
        with open(self._read_processing_stats_path) as read_processing_stats_fh:
            self._read_processing_stats = json.loads(
                read_processing_stats_fh.read()
            )
        with open(self._read_aligner_stats_path) as read_aligner_stats_fh:
            self._read_aligner_stats = json.loads(read_aligner_stats_fh.read())

    def plot_input_read_length(self, plot_path):
        self._plot_read_lengths(
            "read_length_before_processing_and_freq",
            "Length distribution of the input reads - %s",
            plot_path,
        )

    def plot_processed_read_length(self, plot_path):
        self._plot_read_lengths(
            "read_length_after_processing_and_freq",
            "Length distribution of the processed reads - %s",
            plot_path,
        )

    def plot_total_number_of_aligned_reads(self, viz_align_folders_by_species, viz_align_all_folder):
        read_alignment_stats = self._read_in_read_alignment_stats_table()
        all_stats = read_alignment_stats[read_alignment_stats["Species"] == "all"]
        all_stats_selection = all_stats[all_stats["Statistic"].isin(["Total no. of cross aligned reads",
                                                                     "Total no. of unaligned reads"])]
        all_stats_selection.set_index("Statistic", inplace=True)
        all_stats_selection = all_stats_selection.drop(columns=["Species"])
        all_libs = []
        for lib in self._lib_names:
            read_alignment_stats_single_lib = read_alignment_stats[["Species", "Statistic", lib]]
            all_species_exclusive_aligned_reads = []
            for sp, sp_viz_align_folder in viz_align_folders_by_species.items():
                sp_display_name = self._species_folder_prefixes_and_display_names[sp]
                sp_stats = read_alignment_stats_single_lib[read_alignment_stats_single_lib["Species"] == sp_display_name]
                species_exclusive_aligned_reads = sp_stats[sp_stats["Statistic"].isin(["Total no. of uniquely aligned reads",
                                                           "Total no. of split aligned reads",
                                                            "Total no. of multiple aligned reads"])]
                cols = species_exclusive_aligned_reads.columns.to_list()
                exclusive_reads = [sp_display_name, "Total no. of species exclusively aligned reads", species_exclusive_aligned_reads[lib].sum()]
                exclusive_reads_df = pd.DataFrame([exclusive_reads], columns=cols)
                species_exclusive_aligned_reads = pd.concat([species_exclusive_aligned_reads, exclusive_reads_df])
                species_exclusive_aligned_reads.set_index(["Species", "Statistic"], inplace=True)
                all_species_exclusive_aligned_reads.append(species_exclusive_aligned_reads)

            all_species = pd.concat(all_species_exclusive_aligned_reads, axis=0)
            all_libs.append(all_species)
        all_libs_all_species = pd.concat(all_libs, axis=1)
        total_no_of_species_exclusive_reads = all_libs_all_species[all_libs_all_species.index.isin(['Total no. of species exclusively aligned reads'], level=1)]
        total_no_of_species_exclusive_reads.reset_index(inplace=True)

        species_and_statistics_column = total_no_of_species_exclusive_reads["Species"] + " - " + total_no_of_species_exclusive_reads["Statistic"].copy()
        total_no_of_species_exclusive_reads = pd.concat([total_no_of_species_exclusive_reads, species_and_statistics_column], axis=1)
        total_no_of_species_exclusive_reads.drop(["Species", "Statistic"], axis=1, inplace=True)
        total_no_of_species_exclusive_reads.rename(columns={0: "Statistic"}, inplace=True)
        total_no_of_species_exclusive_reads.set_index("Statistic", inplace=True)

        species_and_cross_aligned_and_unaligned = pd.concat([all_stats_selection, total_no_of_species_exclusive_reads])
        print(species_and_cross_aligned_and_unaligned)
        species_and_cross_aligned_and_unaligned_transposed = species_and_cross_aligned_and_unaligned.transpose()
        print(species_and_cross_aligned_and_unaligned_transposed)
        print(viz_align_all_folder)
        sns.set()
        ax = species_and_cross_aligned_and_unaligned_transposed.plot(kind="bar", stacked=True)
        fig = ax.get_figure()
        fig.savefig(f"{viz_align_all_folder}/stacked_species.pdf", bbox_inches='tight')

    def _read_in_read_alignment_stats_table(self):
        print(self._read_alignment_stats_table_path)
        read_alignment_stats = pd.read_csv(
                self._read_alignment_stats_table_path, sep="\t"
            )
        return read_alignment_stats

    def _plot_read_lengths(self, dict_key, title_template, output_file):
        pp = PdfPages(output_file)
        for lib in self._lib_names:
            lengths_and_freqs = self._read_processing_stats[lib][dict_key]
            fig = self._generate_histogram(
                lengths_and_freqs, title_template % lib
            )
            pp.savefig()
            plt.close(fig)
        pp.close()

    def _generate_histogram(self, lengths_and_freqs, title):
        fig = plt.figure()
        lengths = np.array([int(length) for length in lengths_and_freqs.keys()])
        freqs = np.array([int(freq) for freq in lengths_and_freqs.values()])
        ax = fig.add_subplot(111)
        plt.title(title)
        plt.xlabel("Read length [nt]")
        plt.ylabel("Frequency")
        font = {"size": 8}
        matplotlib.rc("font", **font)
        ax.xaxis.set_ticks_position("bottom")
        plt.xticks(np.arange(0, max(lengths) + 1, 10.0))
        plt.bar(lengths, freqs, align="center", color="black", edgecolor="none")
        plt.xlim([0, max(lengths) + 1])
        return fig
