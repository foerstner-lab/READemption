from collections import defaultdict
import json
import numpy as np
import matplotlib
import pandas as pd
import seaborn as sns

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
        viz_align_aligned_reads_by_species_paths,
        species_folder_prefixes_and_display_names,
    ):
        self._lib_names = lib_names
        self._read_processing_stats_path = read_processing_stats_path
        self._read_aligner_stats_path = read_aligner_stats_path
        self._read_alignment_stats_table_path = read_alignment_stats_table_path
        self._viz_align_aligned_reads_by_species_paths = (
            viz_align_aligned_reads_by_species_paths
        )
        self._species_folder_prefixes_and_display_names = (
            species_folder_prefixes_and_display_names
        )

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

    def plot_total_number_of_aligned_reads(
        self, viz_align_all_folder: str
    ) -> None:
        """
        Creates a stacked bar plot of the alignment stats by species and saves
        it as a pdf.
        :param viz_align_all_folder: path of the output folder where the
                pdf for all species will be saved
        """
        # read in the alignment stats table
        read_alignment_stats = self._read_in_read_alignment_stats_table()
        # get statistics for all species combined and
        # select 'Total no. of cross aligned reads' and
        # 'Total no. of unaligned reads'
        all_stats_selection = self._get_all_stats_selection(
            read_alignment_stats
        )
        # get the sums of all species exclusively aligned reads in
        # a combined dataframe.
        total_no_of_species_exclusive_reads = (
            self._get_species_exclusive_reads_combined(read_alignment_stats)
        )

        # Combine the species exclusive stats and the combined stats for
        # all species:
        #                                                     library_one  library_two
        # Statistic
        # Total no. of unaligned reads                                1.0          2.0
        # Total no. of cross aligned reads                            3.0          6.0
        # homo sapiens - Total no. of species exclusively...          7.0         14.0
        # Staphylococcus aureus - Total no. of species ex...          2.0          4.0
        # Influenza A - Total no. of species exclusively ...          1.0          2.0

        species_and_cross_aligned_and_unaligned = pd.concat(
            [all_stats_selection, total_no_of_species_exclusive_reads]
        )
        # Transpose the dataframe to make each lib a bar
        species_and_cross_aligned_and_unaligned_transposed = (
            species_and_cross_aligned_and_unaligned.transpose()
        )
        self._plot_combined_species_align_stats(
            species_and_cross_aligned_and_unaligned_transposed,
            f"{viz_align_all_folder}/stacked_species",
        )
        self._write_combined_species_align_stats(
            species_and_cross_aligned_and_unaligned_transposed,
            f"{viz_align_all_folder}/stacked_species",
        )

    def _write_combined_species_align_stats(
        self, combined_species_aligned_reads_transposed, output_path
    ) -> None:
        combined_species_aligned_reads_transposed.to_csv(
            f"{output_path}.csv", sep="\t"
        )

    def _get_all_stats_selection(
        self, read_alignment_stats: pd.core.frame.DataFrame
    ) -> pd.core.frame.DataFrame:
        """
        Make a pandas Dataframe with the "Total no. of cross aligend reads" and
        the "Total no. of unaligned reads" for the combined species
        :param read_alignment_stats:
        :return: DataFrame e.g.:
                                                      library_one  library_two
            Statistic
            Total no. of unaligned reads              1.0          2.0
            Total no. of cross aligned reads          3.0          6.0
        """
        all_stats = read_alignment_stats[
            read_alignment_stats["Species"] == "all"
        ]
        all_stats_selection = all_stats[
            all_stats["Statistic"].isin(
                [
                    "Total no. of cross aligned reads",
                    "Total no. of unaligned reads",
                ]
            )
        ]
        # move statistic column to index and remove species information:
        #                                              library_one  library_two
        #    Statistic
        #    Total no. of unaligned reads              1.0          2.0
        #    Total no. of cross aligned reads          3.0          6.0

        all_stats_selection.set_index("Statistic", inplace=True)
        all_stats_selection = all_stats_selection.drop(columns=["Species"])
        return all_stats_selection

    def _get_species_exclusive_reads_combined(
        self, read_alignment_stats: pd.core.frame.DataFrame
    ) -> pd.core.frame.DataFrame:
        """
        create and return a dataframe that contains the sums of the species
        exclusively aligned reads of each species
        :param read_alignment_stats:
        :return: total_no_of_species_exclusive_reads:
        """
        all_species_by_lib = []
        for lib in self._lib_names:
            # get the alignment stats for a single lib
            read_alignment_stats_single_lib = read_alignment_stats[
                ["Species", "Statistic", lib]
            ]
            all_species_exclusive_aligned_reads = []
            # get the reads that are exlusively aligned to each species
            for (
                sp,
                sp_display_name,
            ) in self._species_folder_prefixes_and_display_names.items():
                # select statistics by species display name
                sp_stats = read_alignment_stats_single_lib[
                    read_alignment_stats_single_lib["Species"]
                    == sp_display_name
                ]
                ## Create Dataframe for the combined species
                # select the statistics that only contain exclusively aligned
                # reads for the given species
                species_exclusive_aligned_reads = sp_stats[
                    sp_stats["Statistic"].isin(
                        [
                            "Total no. of uniquely aligned reads",
                            "Total no. of split aligned reads",
                            "Total no. of multiple aligned reads",
                        ]
                    )
                ]
                # retrieve the column names e.g.:
                # cols = ['Species', 'Statistic', 'library_one']
                cols = species_exclusive_aligned_reads.columns.to_list()
                # get a list that contains the sum of the species exclusively
                # aligned reads e.g.:
                # exculisve_reads =
                #   ['homo sapiens',
                #   'Total no. of species exclusively aligned reads',
                #    7.0]
                exclusive_reads = [
                    sp_display_name,
                    "Total no. of species exclusively aligned reads",
                    species_exclusive_aligned_reads[lib].sum(),
                ]
                # make dataframe of species exclusive reads e.g.:
                #        Species      Statistic                                          library_two
                #   0   Influenza A   Total no. of species exclusively aligned reads     2.0
                exclusive_reads_df = pd.DataFrame(
                    [exclusive_reads], columns=cols
                )
                # Combine the species exclusive reads and the sum of them in a
                # dataframe and set 'Species' and 'Statistic' columns as index
                #                                                      library_two
                # Species     Statistic
                # homo sapiens Total no. of uniquely aligned reads             6.0
                #             Total no. of split aligned reads                6.0
                #             Total no. of multiple aligned reads             2.0
                #             Total no. of species exclusively aligned reads  14.0
                species_exclusive_aligned_reads = pd.concat(
                    [species_exclusive_aligned_reads, exclusive_reads_df]
                )

                species_exclusive_aligned_reads.set_index(
                    ["Species", "Statistic"], inplace=True
                )
                # store the dataframes of all species in a list
                all_species_exclusive_aligned_reads.append(
                    species_exclusive_aligned_reads
                )
            # combine all dataframes of all species of a single lib
            all_species = pd.concat(all_species_exclusive_aligned_reads, axis=0)
            # safe all dataframes of a lib to a list
            all_species_by_lib.append(all_species)

        ## Create Dataframe for the combined species
        # combine all dataframes of all libs (all species)
        all_libs_all_species = pd.concat(all_species_by_lib, axis=1)

        # select only the sum of the species exclusively aligned reads
        total_no_of_species_exclusive_reads = all_libs_all_species[
            all_libs_all_species.index.isin(
                ["Total no. of species exclusively aligned reads"], level=1
            )
        ]
        # reset index and return species to columns:
        #                 Species                                       Statistic  library_one  library_two
        # 0           homo sapiens  Total no. of species exclusively aligned reads          7.0         14.0
        # 1  Staphylococcus aureus  Total no. of species exclusively aligned reads          2.0          4.0
        # 2            Influenza A  Total no. of species exclusively aligned reads          1.0          2.0
        total_no_of_species_exclusive_reads.reset_index(inplace=True)
        # Make new Series with 'Species' and 'Statistic' connected with ' - ':
        # 0    homo sapiens - Total no. of species exclusivel...
        # 1    Staphylococcus aureus - Total no. of species e...
        # 2    Influenza A - Total no. of species exclusively...

        species_and_statistics_column = (
            total_no_of_species_exclusive_reads["Species"]
            + " - "
            + total_no_of_species_exclusive_reads["Statistic"].copy()
        )

        #####
        # total_no_of_species_exclusive_reads["Species Statistic"] = total_no_of_species_exclusive_reads["Species"] + " - " + total_no_of_species_exclusive_reads["Statistic"]
        # print(total_no_of_species_exclusive_reads["Species Statistic"])
        # print(total_no_of_species_exclusive_reads)
        #####

        # Add the new Series as a column
        #                  Species                                       Statistic  library_one  library_two                                                  0
        # 0           homo sapiens  Total no. of species exclusively aligned reads          7.0         14.0  homo sapiens - Total no. of species exclusivel...
        # 1  Staphylococcus aureus  Total no. of species exclusively aligned reads          2.0          4.0  Staphylococcus aureus - Total no. of species e...
        # 2            Influenza A  Total no. of species exclusively aligned reads          1.0          2.0  Influenza A - Total no. of species exclusively...

        total_no_of_species_exclusive_reads = pd.concat(
            [
                total_no_of_species_exclusive_reads,
                species_and_statistics_column,
            ],
            axis=1,
        )
        # Drop the 'Species' and 'Statistic' columns and rename the combined
        # column to 'Statistic'. Then set it as index:
        #                                                     library_one  library_two
        # Statistic
        # homo sapiens - Total no. of species exclusively...          7.0         14.0
        # Staphylococcus aureus - Total no. of species ex...          2.0          4.0
        # Influenza A - Total no. of species exclusively ...          1.0          2.0

        total_no_of_species_exclusive_reads.drop(
            ["Species", "Statistic"], axis=1, inplace=True
        )
        total_no_of_species_exclusive_reads.rename(
            columns={0: "Statistic"}, inplace=True
        )
        total_no_of_species_exclusive_reads.set_index("Statistic", inplace=True)
        return total_no_of_species_exclusive_reads

    def _plot_combined_species_align_stats(
        self,
        species_and_cross_aligned_and_unaligned_transposed: pd.core.frame.DataFrame,
        output_path: str,
    ) -> None:
        """
        :param species_and_cross_aligned_and_unaligned_transposed: Dataframe
        that contains the species aligned, cross aligned and unaligned reads
        for all species
        :param output_path: the path where to save the figure
        """
        sns.set_theme(palette="colorblind")
        ax = species_and_cross_aligned_and_unaligned_transposed.plot(
            kind="bar", stacked=True
        )
        # reverse the labels to have them in the same order as the bar segments
        handles, labels = plt.gca().get_legend_handles_labels()
        handles.reverse()
        labels.reverse()
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.35) ,handles=handles, labels=labels)

        fig = ax.get_figure()
        fig.savefig(f"{output_path}.pdf", bbox_inches="tight")

    def plot_species_exclusive_reads_for_each_species(self) -> None:
        """
        Creates bar plots with different alignment stats for
        each species separately
        """
        # read in the alignment stats table
        read_alignment_stats = self._read_in_read_alignment_stats_table()
        ## Create Dataframes for single species
        aligned_reads_by_species = defaultdict(list)
        for lib in self._lib_names:
            # get the alignment stats for a single lib
            read_alignment_stats_single_lib = read_alignment_stats[
                ["Species", "Statistic", lib]
            ]
            all_species_exclusive_aligned_reads = []
            # get the reads that are exlusively aligned to each species
            for (
                sp,
                sp_display_name,
            ) in self._species_folder_prefixes_and_display_names.items():
                # select statistics by species display name
                sp_stats = read_alignment_stats_single_lib[
                    read_alignment_stats_single_lib["Species"]
                    == sp_display_name
                ]

                ## Create Dataframes for single species
                # select all kind of aligned reads of a species
                species_aligned_reads = sp_stats[
                    sp_stats["Statistic"].isin(
                        [
                            "Total no. of uniquely aligned reads",
                            "Total no. of split aligned reads",
                            "Total no. of multiple aligned reads",
                            "Total no. of cross aligned reads",
                        ]
                    )
                ]
                # Drop the 'Species' column
                species_aligned_reads = species_aligned_reads.drop(
                    "Species", axis=1
                )
                # Move the 'Statistic' column to index
                species_aligned_reads.set_index("Statistic", inplace=True)
                # Add the library stats to the species stats
                aligned_reads_by_species[sp].append(species_aligned_reads)
        for (
            sp,
            lib_stats,
        ) in aligned_reads_by_species.items():
            species_aligned_reads = pd.concat(lib_stats, axis=1)
            species_aligned_reads_transposed = species_aligned_reads.transpose()
            output_path = self._viz_align_aligned_reads_by_species_paths[sp]
            self._plot_species_align_stats(
                species_aligned_reads_transposed, output_path
            )
            self._write_species_align_stats(
                species_aligned_reads_transposed, output_path
            )

    def _write_species_align_stats(
        self, species_aligned_reads_transposed, output_path
    ) -> None:
        species_aligned_reads_transposed.to_csv(f"{output_path}.csv", sep="\t")

    def _plot_species_align_stats(
        self, species_aligned_reads_transposed, output_path
    ) -> None:
        sns.set_theme(palette="colorblind")
        ax = species_aligned_reads_transposed.plot(kind="bar", stacked=True)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.35))
        fig = ax.get_figure()
        fig.savefig(f"{output_path}.pdf", bbox_inches="tight")

    def _read_in_read_alignment_stats_table(self):
        return pd.read_csv(self._read_alignment_stats_table_path, sep="\t")

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
