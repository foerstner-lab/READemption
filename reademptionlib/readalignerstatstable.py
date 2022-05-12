import pandas as pd


class ReadAlignerStatsTable(object):
    def __init__(
        self,
        read_processing_stats,
        alignment_stats,
        libs,
        output_path_table,
        output_path_table_transposed,
        paired_end,
        species_folder_prefixes_and_display_names,
        references_by_species,
        fragments=False
    ):
        self._table = []
        self._read_processing_stats = read_processing_stats
        self._alignment_stats = alignment_stats
        self._libs = libs
        self._output_path_table = output_path_table
        self._output_path_table_transposed = output_path_table_transposed
        self._paired_end = paired_end
        self._species_folder_prefixes_and_display_names = (
            species_folder_prefixes_and_display_names
        )
        self._references_by_species = references_by_species
        self._fragments = fragments

        if self._fragments:
            self.reads_or_fragments = "fragments"
        else:
            self.reads_or_fragments = "reads"

    def write(self):
        all_stats = self._create_table_all_statistics()
        all_stats.to_csv(self._output_path_table, sep="\t")
        all_stats.transpose().to_csv(
            self._output_path_table_transposed, sep="\t"
        )

    def _create_table_all_statistics(self):
        total_stats = self._create_statistics_table_total()
        species_stats = self._create_statistics_table_species()
        reference_stats = self._reference_wise_stats()
        all_stats = pd.concat([total_stats, species_stats, reference_stats])
        return all_stats

    def _create_statistics_table_total(self):
        stats_total_all = []
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
                        f"no_of_aligned_{self.reads_or_fragments}"
                    ]
                ),
                f"Total no. of aligned {self.reads_or_fragments}",
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
                        f"no_of_uniquely_aligned_{self.reads_or_fragments}"
                    ]
                ),
                f"Total no. of uniquely aligned {self.reads_or_fragments}",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                int(
                    self._alignment_stats[lib]["stats_total"][
                        f"no_of_split_aligned_{self.reads_or_fragments}"
                    ]
                ),
                f"Total no. of split aligned {self.reads_or_fragments}",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                int(
                    self._alignment_stats[lib]["stats_total"][
                        f"no_of_multiple_aligned_{self.reads_or_fragments}"
                    ]
                ),
                f"Total no. of multiple aligned {self.reads_or_fragments}",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                int(
                    self._alignment_stats[lib]["stats_total"][
                        f"no_of_cross_aligned_{self.reads_or_fragments}"
                    ]
                ),
                f"Total no. of cross aligned {self.reads_or_fragments}",
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
            # Do not add to fragment stats
            if not self._fragments:
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
            # Do not add to fragment stats
            if not self._fragments:
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
                                f"no_of_uniquely_aligned_{self.reads_or_fragments}"
                            ]
                        ),
                        (
                            self._alignment_stats[lib]["stats_total"][
                                f"no_of_aligned_{self.reads_or_fragments}"
                            ]
                        ),
                    ),
                    2,
                ),
                f"Percentage of uniquely aligned {self.reads_or_fragments} (in relation to all aligned {self.reads_or_fragments})",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                round(
                    self._calc_percentage(
                        (
                            self._alignment_stats[lib]["stats_total"][
                                f"no_of_split_aligned_{self.reads_or_fragments}"
                            ]
                        ),
                        (
                            self._alignment_stats[lib]["stats_total"][
                                f"no_of_aligned_{self.reads_or_fragments}"
                            ]
                        ),
                    ),
                    2,
                ),
                f"Percentage of split aligned {self.reads_or_fragments} (in relation to all aligned {self.reads_or_fragments})",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                round(
                    self._calc_percentage(
                        (
                            self._alignment_stats[lib]["stats_total"][
                                f"no_of_multiple_aligned_{self.reads_or_fragments}"
                            ]
                        ),
                        (
                            self._alignment_stats[lib]["stats_total"][
                                f"no_of_aligned_{self.reads_or_fragments}"
                            ]
                        ),
                    ),
                    2,
                ),
                f"Percentage of multiple aligned {self.reads_or_fragments} (in relation to all aligned {self.reads_or_fragments})",
            )
            stats_total = self._append_to_df(
                stats_total,
                lib,
                round(
                    self._calc_percentage(
                        (
                            self._alignment_stats[lib]["stats_total"][
                                f"no_of_cross_aligned_{self.reads_or_fragments}"
                            ]
                        ),
                        (
                            self._alignment_stats[lib]["stats_total"][
                                f"no_of_aligned_{self.reads_or_fragments}"
                            ]
                        ),
                    ),
                    2,
                ),
                f"Percentage of cross aligned {self.reads_or_fragments} (in relation to all aligned {self.reads_or_fragments})",
            )
            stats_total.insert(0, "Species", "all")
            stats_total["Statistic"] = stats_total.index
            stats_total.reset_index(drop=True, inplace=True)
            stats_total.set_index(["Species", "Statistic"], inplace=True)
            stats_total_all.append(stats_total)
        stats_total_combined = pd.concat(stats_total_all, axis=1)
        return stats_total_combined

    def _create_statistics_table_species(self):
        species_tables = []
        for (
            sp,
            sp_display_name,
        ) in self._species_folder_prefixes_and_display_names.items():
            stats_all_libs = []
            for lib in self._libs:
                lib_stats = self._create_overview_stats(
                    lib, sp, sp_display_name
                )
                lib_stats.set_index(["Species", "Statistic"], inplace=True)
                stats_all_libs.append(lib_stats)
            sp_stats = pd.concat(stats_all_libs, axis=1)
            species_tables.append(sp_stats)
        combined_species_table = pd.concat(species_tables)
        return combined_species_table

    def _create_overview_stats(self, lib, species, species_display_name):
        stats_total = pd.DataFrame(columns=[lib])
        stats_total = self._append_to_df(
            stats_total,
            lib,
            int(
                self._alignment_stats[lib]["species_stats"][species][
                    f"no_of_aligned_{self.reads_or_fragments}"
                ]
            ),
            f"Total no. of aligned {self.reads_or_fragments}",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            int(
                self._alignment_stats[lib]["species_stats"][species][
                    f"no_of_uniquely_aligned_{self.reads_or_fragments}"
                ]
            ),
            f"Total no. of uniquely aligned {self.reads_or_fragments}",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            int(
                self._alignment_stats[lib]["species_stats"][species][
                    f"no_of_split_aligned_{self.reads_or_fragments}"
                ]
            ),
            f"Total no. of split aligned {self.reads_or_fragments}",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            int(
                self._alignment_stats[lib]["species_stats"][species][
                    f"no_of_multiple_aligned_{self.reads_or_fragments}"
                ]
            ),
            f"Total no. of multiple aligned {self.reads_or_fragments}",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            int(
                self._alignment_stats[lib]["species_stats"][species][
                    f"no_of_cross_aligned_{self.reads_or_fragments}"
                ]
            ),
            f"Total no. of cross aligned {self.reads_or_fragments}",
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
        # Do not add to fragment stats
        if not self._fragments:
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
        # Do not add to fragment stats
        if not self._fragments:
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
                            f"no_of_uniquely_aligned_{self.reads_or_fragments}"
                        ]
                    ),
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            f"no_of_aligned_{self.reads_or_fragments}"
                        ]
                    ),
                ),
                2,
            ),
            f"Percentage of uniquely aligned {self.reads_or_fragments} (in relation to all aligned {self.reads_or_fragments})",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            round(
                self._calc_percentage(
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            f"no_of_split_aligned_{self.reads_or_fragments}"
                        ]
                    ),
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            f"no_of_aligned_{self.reads_or_fragments}"
                        ]
                    ),
                ),
                2,
            ),
            f"Percentage of split aligned {self.reads_or_fragments} (in relation to all aligned {self.reads_or_fragments})",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            round(
                self._calc_percentage(
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            f"no_of_multiple_aligned_{self.reads_or_fragments}"
                        ]
                    ),
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            f"no_of_aligned_{self.reads_or_fragments}"
                        ]
                    ),
                ),
                2,
            ),
            f"Percentage of multiple aligned {self.reads_or_fragments} (in relation to all aligned {self.reads_or_fragments})",
        )
        stats_total = self._append_to_df(
            stats_total,
            lib,
            round(
                self._calc_percentage(
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            f"no_of_cross_aligned_{self.reads_or_fragments}"
                        ]
                    ),
                    (
                        self._alignment_stats[lib]["species_stats"][species][
                            f"no_of_aligned_{self.reads_or_fragments}"
                        ]
                    ),
                ),
                2,
            ),
            f"Percentage of cross aligned {self.reads_or_fragments} (in relation to all aligned {self.reads_or_fragments})",
        )
        stats_total.insert(0, "Species", species_display_name)
        stats_total["Statistic"] = stats_total.index
        stats_total.reset_index(drop=True, inplace=True)
        stats_total = stats_total[["Species", "Statistic", lib]]
        return stats_total

    def _reference_wise_stats(self):
        stats_per_ref_all_libs = []
        for lib in self._libs:
            lib_data_frames = []
            for species, ref_ids in self._references_by_species.items():
                species_for_ref_id = species
                for ref_id in ref_ids:

                    species_display_name = (
                        self._species_folder_prefixes_and_display_names[
                            species_for_ref_id
                        ]
                    )
                    ref_stats = pd.DataFrame(columns=[lib])
                    ref_stats = self._append_to_df(
                        ref_stats,
                        lib,
                        int(
                            self._alignment_stats[lib]["stats_per_reference"][
                                ref_id
                            ][f"no_of_aligned_{self.reads_or_fragments}"]
                        ),
                        f"{ref_id} - No. of aligned {self.reads_or_fragments}",
                    )
                    ref_stats = self._append_to_df(
                        ref_stats,
                        lib,
                        int(
                            self._alignment_stats[lib]["stats_per_reference"][
                                ref_id
                            ][f"no_of_uniquely_aligned_{self.reads_or_fragments}"]
                        ),
                        f"{ref_id} - No. of uniquely aligned {self.reads_or_fragments}",
                    )
                    ref_stats = self._append_to_df(
                        ref_stats,
                        lib,
                        int(
                            self._alignment_stats[lib]["stats_per_reference"][
                                ref_id
                            ][f"no_of_split_aligned_{self.reads_or_fragments}"]
                        ),
                        f"{ref_id} - No. of split aligned {self.reads_or_fragments}",
                    )
                    ref_stats = self._append_to_df(
                        ref_stats,
                        lib,
                        int(
                            self._alignment_stats[lib]["stats_per_reference"][
                                ref_id
                            ][f"no_of_multiple_aligned_{self.reads_or_fragments}"]
                        ),
                        f"{ref_id} - No. of multiple aligned {self.reads_or_fragments}",
                    )
                    ref_stats = self._append_to_df(
                        ref_stats,
                        lib,
                        int(
                            self._alignment_stats[lib]["stats_per_reference"][
                                ref_id
                            ]["no_of_alignments"]
                        ),
                        f"{ref_id} - No. of alignments",
                    )
                    ref_stats.insert(0, "Species", species_display_name)
                    ref_stats.reset_index(inplace=True)
                    ref_stats.rename(
                        columns={"index": "Statistic"}, inplace=True
                    )
                    ref_stats.set_index(["Species", "Statistic"], inplace=True)

                    lib_data_frames.append(ref_stats)
            stats_per_ref_lib = pd.concat(lib_data_frames)
            stats_per_ref_all_libs.append(stats_per_ref_lib)
        stats_per_ref_combined = pd.concat(stats_per_ref_all_libs, axis=1)
        return stats_per_ref_combined

    def _append_to_df(self, df, lib, stat, value):
        dataframe_to_add = pd.DataFrame({lib: stat}, index=[value])
        df = pd.concat([df, dataframe_to_add], axis=0)
        return df

    def _get_read_process_number(self, lib, attribute):
        factor = 1
        if self._paired_end:
            factor = 2
        return self._read_processing_stats[lib][attribute] * factor

    def _calc_percentage(self, mult, div):
        try:
            return float(mult) / float(div) * 100
        except ZeroDivisionError:
            return 0.0
