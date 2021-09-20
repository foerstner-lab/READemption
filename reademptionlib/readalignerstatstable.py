import pandas as pd
import pprint


class ReadAlignerStatsTable(object):
    def __init__(
        self,
        read_processing_stats,
        alignment_stats,
        libs,
        output_path,
        paired_end,
        species_folder_prefixes_and_display_names
    ):
        self._table = []
        self._read_processing_stats = read_processing_stats
        self._alignment_stats = alignment_stats
        self._libs = libs
        self._output_path = output_path
        self._paired_end = paired_end
        self._species_folder_prefixes_and_display_names = species_folder_prefixes_and_display_names


        self._create_statistics_table_total()
        self._create_statistics_table_species()

    def _create_statistics_table_total(self):
        # print(self._alignment_stats)
        pprint.pprint(self._read_processing_stats)
        for lib in self._libs:
            stats_total = pd.DataFrame(columns=[lib])
            stats_total = self._append_to_df(
                stats_total,
                lib,
                self._get_read_process_number(lib, "total_no_of_reads"),
                "No. of input reads",
            )

            stats_total = self._append_to_df(
                stats_total,
                lib,
                self._get_read_process_number(lib, "polya_removed"),
                "No. of reads - PolyA detected and removed",
            )

            stats_total = self._append_to_df(
                stats_total,
                lib,
                self._get_read_process_number(lib, "single_a_removed"),
                "No. of reads - Single 3' A removed",
            )

            stats_total = self._append_to_df(
                stats_total,
                lib,
                self._get_read_process_number(lib, "unmodified"),
                "No. of reads - Unmodified",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                self._get_read_process_number(lib, "too_short"),
                "No. of reads - Removed as too short",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                self._get_read_process_number(lib, "long_enough"),
                "No. of reads - Long enough and used for alignment",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                int(
                    self._alignment_stats[lib]["stats_total"][
                        "no_of_aligned_reads"
                    ]
                ),
                "Total no. of aligned reads",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                int(
                    self._alignment_stats[lib]["stats_total"][
                        "no_of_unaligned_reads"
                    ]
                ),
                "Total no. of unaligned reads",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                int(
                    self._alignment_stats[lib]["stats_total"][
                        "no_of_uniquely_aligned_reads"
                    ]
                ),
                "Total no. of uniquely aligned reads",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                int(
                    self._alignment_stats[lib]["stats_total"][
                        "no_of_split_aligned_reads"
                    ]
                ),
                "Total no. of split aligned reads",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                int(
                    self._alignment_stats[lib]["stats_total"][
                        "no_of_multiple_aligned_reads"
                    ]
                ),
                "Total no. of multiple aligned reads",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                int(
                    self._alignment_stats[lib]["stats_total"][
                        "no_of_cross_aligned_reads"
                    ]
                ),
                "Total no. of cross aligned reads",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                int(
                    self._alignment_stats[lib]["stats_total"][
                        "no_of_alignments"
                    ]
                ),
                "Total no. of alignments",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                round(
                    self._calc_percentage(
                        (
                            self._alignment_stats[lib]["stats_total"][
                                "no_of_aligned_reads"
                            ]
                        ),
                        self._get_read_process_number(lib, "total_no_of_reads"),
                    ),
                    2,
                ),
                "Percentage of aligned reads (compared to no. of input reads)",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                round(
                    self._calc_percentage(
                        (
                            self._alignment_stats[lib]["stats_total"][
                                "no_of_aligned_reads"
                            ]
                        ),
                        self._get_read_process_number(lib, "long_enough"),
                    ),
                    2,
                ),
                "Percentage of aligned reads (compared to no. of long enough reads)",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                round(
                    self._calc_percentage(
                        (
                            self._alignment_stats[lib]["stats_total"][
                                "no_of_uniquely_aligned_reads"
                            ]
                        ),
                        (
                            self._alignment_stats[lib]["stats_total"][
                                "no_of_aligned_reads"
                            ]
                        ),
                    ),
                    2,
                ),
                "Percentage of uniquely aligned reads (in relation to all aligned reads)",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                round(
                    self._calc_percentage(
                        (
                            self._alignment_stats[lib]["stats_total"][
                                "no_of_split_aligned_reads"
                            ]
                        ),
                        (
                            self._alignment_stats[lib]["stats_total"][
                                "no_of_aligned_reads"
                            ]
                        ),
                    ),
                    2,
                ),
                "Percentage of split aligned reads (in relation to all aligned reads)",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                round(
                    self._calc_percentage(
                        (
                            self._alignment_stats[lib]["stats_total"][
                                "no_of_multiple_aligned_reads"
                            ]
                        ),
                        (
                            self._alignment_stats[lib]["stats_total"][
                                "no_of_aligned_reads"
                            ]
                        ),
                    ),
                    2,
                ),
                "Percentage of multiple aligned reads (in relation to all aligned reads)",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                round(
                    self._calc_percentage(
                        (
                            self._alignment_stats[lib]["stats_total"][
                                "no_of_cross_aligned_reads"
                            ]
                        ),
                        (
                            self._alignment_stats[lib]["stats_total"][
                                "no_of_aligned_reads"
                            ]
                        ),
                    ),
                    2,
                ),
                "Percentage of cross aligned reads (in relation to all aligned reads)",
            )
            stats_total.insert(0, 'Species', "all")
            stats_total['Statistic'] = stats_total.index
            stats_total.reset_index(drop=True, inplace=True)
            stats_total = stats_total[["Species", "Statistic", lib]]

        #print(stats_total)

    def _create_statistics_table_species(self):
        species_tables = []
        for sp, sp_display_name in self._species_folder_prefixes_and_display_names.items():
            for lib in self._libs:
                stats = self._create_overview_stats(lib, sp, sp_display_name)
                species_tables.append(stats)
        combined_species_table = pd.concat(species_tables)
        print(combined_species_table)


    def _create_overview_stats(self, lib, species, species_display_name):
        stats_total = pd.DataFrame(columns=[lib])
        stats_total = self._append_to_df(
            stats_total,
            lib,
            int(
                self._alignment_stats[lib]["species_stats"][species][
                    "no_of_aligned_reads"
                ]
            ),
            "Total no. of aligned reads",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            int(
                self._alignment_stats[lib]["species_stats"][species][
                    "no_of_uniquely_aligned_reads"
                ]
            ),
            "Total no. of uniquely aligned reads",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            int(
                self._alignment_stats[lib]["species_stats"][species][
                    "no_of_split_aligned_reads"
                ]
            ),
            "Total no. of split aligned reads",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            int(
                self._alignment_stats[lib]["species_stats"][species][
                    "no_of_multiple_aligned_reads"
                ]
            ),
            "Total no. of multiple aligned reads",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            int(
                self._alignment_stats[lib]["species_stats"][species][
                    "no_of_cross_aligned_reads"
                ]
            ),
            "Total no. of cross aligned reads",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            int(
                self._alignment_stats[lib]["species_stats"][species][
                    "no_of_alignments"
                ]
            ),
            "Total no. of alignments",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            round(
                self._calc_percentage(
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            "no_of_aligned_reads"
                        ]
                    ),
                    self._get_read_process_number(lib, "total_no_of_reads"),
                ),
                2,
            ),
            "Percentage of aligned reads (compared to no. of input reads)",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            round(
                self._calc_percentage(
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            "no_of_aligned_reads"
                        ]
                    ),
                    self._get_read_process_number(lib, "long_enough"),
                ),
                2,
            ),
            "Percentage of aligned reads (compared to no. of long enough reads)",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            round(
                self._calc_percentage(
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            "no_of_uniquely_aligned_reads"
                        ]
                    ),
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            "no_of_aligned_reads"
                        ]
                    ),
                ),
                2,
            ),
            "Percentage of uniquely aligned reads (in relation to all aligned reads)",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            round(
                self._calc_percentage(
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            "no_of_split_aligned_reads"
                        ]
                    ),
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            "no_of_aligned_reads"
                        ]
                    ),
                ),
                2,
            ),
            "Percentage of split aligned reads (in relation to all aligned reads)",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            round(
                self._calc_percentage(
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            "no_of_multiple_aligned_reads"
                        ]
                    ),
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            "no_of_aligned_reads"
                        ]
                    ),
                ),
                2,
            ),
            "Percentage of multiple aligned reads (in relation to all aligned reads)",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            round(
                self._calc_percentage(
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            "no_of_cross_aligned_reads"
                        ]
                    ),
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            "no_of_aligned_reads"
                        ]
                    ),
                ),
                2,
            ),
            "Percentage of cross aligned reads (in relation to all aligned reads)",
        )
        stats_total.insert(0, 'Species', species_display_name)
        stats_total['Statistic'] = stats_total.index
        stats_total.reset_index(drop=True, inplace=True)
        stats_total = stats_total[["Species", "Statistic", lib]]
        return stats_total


    def _append_to_df(self, df, lib, stat, value):
        df = df.append(pd.Series({lib: stat}, name=value))
        return df

    def write(self):
        self._add_global_countings()
        self._add_reference_wise_coutings()
        with open(self._output_path, "w") as table_fh:
            table_fh.write(
                "\n".join(
                    [
                        "\t".join([str(cell) for cell in row])
                        for row in self._table
                    ]
                )
                + "\n"
            )

    def _add_global_countings(self):
        for title, data in [
            ("Libraries", self._libs),
            (
                "No. of input reads",
                self._get_read_process_numbers("total_no_of_reads"),
            ),
            (
                "No. of reads - PolyA detected and removed",
                self._get_read_process_numbers("polya_removed"),
            ),
            (
                "No. of reads - Single 3' A removed",
                self._get_read_process_numbers("single_a_removed"),
            ),
            (
                "No. of reads - Unmodified",
                self._get_read_process_numbers("unmodified"),
            ),
            (
                "No. of reads - Removed as too short",
                self._get_read_process_numbers("too_short"),
            ),
            (
                "No. of reads - Long enough and used for alignment",
                self._get_read_process_numbers("long_enough"),
            ),
            (
                "Total no. of aligned reads",
                self._total_alignment_stat_numbers("no_of_aligned_reads"),
            ),
            (
                "Total no. of unaligned reads",
                self._total_alignment_stat_numbers("no_of_unaligned_reads"),
            ),
            (
                "Total no. of uniquely aligned reads",
                self._total_alignment_stat_numbers(
                    "no_of_uniquely_aligned_reads"
                ),
            ),
            (
                "Total no. of alignments",
                self._total_alignment_stat_numbers("no_of_alignments"),
            ),
            (
                "Total no. of split aligned reads",
                self._total_alignment_stat_numbers("no_of_split_aligned_reads"),
            ),
            (
                "Total no. of multiple aligned reads",
                self._total_alignment_stat_numbers(
                    "no_of_multiple_aligned_reads"
                ),
            ),
            (
                "Percentage of aligned reads (compared to no. of input reads)",
                self._perc_aligned_reads_all_input(),
            ),
            (
                "Percentage of aligned reads (compared to no. of long enough "
                "reads)",
                self._perc_aligned_reads_all_long_enough(),
            ),
            (
                "Percentage of uniquely aligned reads (in relation to all aligned"
                " reads)",
                self._perc_uniquely_aligned_reads(),
            ),
        ]:
            self._table.append([title] + data)

    def _add_reference_wise_coutings(self):
        ref_ids = sorted(
            list(
                list(self._alignment_stats.values())[0][
                    "stats_per_reference"
                ].keys()
            )
        )
        for ref_id in ref_ids:
            for title_template, data in [
                (
                    "%s - No. of aligned reads",
                    self._alignment_number_per_ref_seq(
                        ref_id, "no_of_aligned_reads"
                    ),
                ),
                (
                    "%s - No. of uniquely aligned reads",
                    self._alignment_number_per_ref_seq(
                        ref_id, "no_of_uniquely_aligned_reads"
                    ),
                ),
                (
                    "%s - No. of split aligned reads",
                    self._alignment_number_per_ref_seq(
                        ref_id, "no_of_split_aligned_reads"
                    ),
                ),
                (
                    "%s - No. of multiple aligned reads",
                    self._alignment_number_per_ref_seq(
                        ref_id, "no_of_multiple_aligned_reads"
                    ),
                ),
                (
                    "%s - No. of alignments",
                    self._alignment_number_per_ref_seq(
                        ref_id, "no_of_alignments"
                    ),
                ),
            ]:
                self._table.append([title_template % ref_id] + data)

    def _alignment_number_per_ref_seq(self, ref_id, attribute):
        return [
            round(
                self._alignment_stats[lib]["stats_per_reference"][ref_id].get(
                    attribute, 0
                )
            )
            for lib in self._libs
        ]

    def _total_alignment_stat_numbers(self, attribute, round_nums=True):
        countings = [
            self._alignment_stats[lib]["stats_total"].get(attribute, 0)
            for lib in self._libs
        ]
        if round_nums is True:
            return [round(counting) for counting in countings]
        else:
            return countings

    # TODO remove _get_read_process_numbers
    def _get_read_process_numbers(self, attribute):
        factor = 1
        if self._paired_end:
            factor = 2
        return [
            self._read_processing_stats[lib][attribute] * factor
            for lib in self._libs
        ]

    def _get_read_process_number(self, lib, attribute):
        factor = 1
        if self._paired_end:
            factor = 2
        return self._read_processing_stats[lib][attribute] * factor

    def _perc_aligned_reads_all_input(self):
        return [
            round(self._calc_percentage(aligned_reads, total_reads), 2)
            for aligned_reads, total_reads in zip(
                self._total_alignment_stat_numbers(
                    "no_of_aligned_reads", round_nums=False
                ),
                self._get_read_process_numbers("total_no_of_reads"),
            )
        ]

    def _perc_aligned_reads_all_long_enough(self):
        return [
            round(self._calc_percentage(aligned_reads, total_reads), 2)
            for aligned_reads, total_reads in zip(
                self._total_alignment_stat_numbers(
                    "no_of_aligned_reads", round_nums=False
                ),
                self._get_read_process_numbers("long_enough"),
            )
        ]

    def _perc_uniquely_aligned_reads(self):
        return [
            round(
                self._calc_percentage(uniquely_aligned_reads, aligned_reads), 2
            )
            for uniquely_aligned_reads, aligned_reads in zip(
                self._total_alignment_stat_numbers(
                    "no_of_uniquely_aligned_reads"
                ),
                self._total_alignment_stat_numbers("no_of_aligned_reads"),
            )
        ]

    def _calc_percentage(self, mult, div):
        try:
            return float(mult) / float(div) * 100
        except ZeroDivisionError:
            return 0.0
