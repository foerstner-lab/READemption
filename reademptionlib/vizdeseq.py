import csv
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd


class DESeqViz(object):
    def __init__(self, deseq_script_path, deseq_path_template, max_pvalue):
        self._deseq_script_path = deseq_script_path
        self._deseq_path_template = deseq_path_template
        self._basemean_column = 1
        self._log_2_fold_chance_column = 2
        self._p_value_column = 5
        self._adj_p_value_column = 6
        self._p_value_significance_limit = max_pvalue
        self._log_fold_change_limit = 2
        self._log_2_fold_chance_limit = 1

    def create_scatter_plots(self, plot_path):
        self._pp_scatterplots = PdfPages(plot_path)
        conditions = self._extract_condition_names()
        for condition_1 in conditions:
            for condition_2 in conditions:
                if condition_1 == condition_2:
                    continue
                self._create_scatter_plots(condition_1, condition_2)
        self._pp_scatterplots.close()

    def _create_scatter_plots(self, condition_1, condition_2):
        matplotlib.style.use("ggplot")
        font = {"family": "sans-serif", "size": 7}
        matplotlib.rc("font", **font)
        deseq_result = pd.read_table(
            self._deseq_path_template + f"{condition_1}_vs_{condition_2}.csv"
        )
        # Remove 0 base mean row as those would otherwise cause trouble
        # for the log10
        deseq_result = deseq_result[deseq_result.baseMean > 0]
        fig = plt.figure()
        significant_deseq_result = deseq_result[
            (deseq_result.padj < self._p_value_significance_limit)
        ]
        non_significant_deseq_result = deseq_result[
            (deseq_result.padj >= self._p_value_significance_limit)
        ]
        plt.plot(
            np.log10(non_significant_deseq_result.baseMean),
            non_significant_deseq_result.log2FoldChange,
            ".k",
            alpha=0.5,
        )
        plt.plot(
            np.log10(significant_deseq_result.baseMean),
            significant_deseq_result.log2FoldChange,
            ".r",
            alpha=0.3,
        )
        plt.title(f"{condition_1} vs. {condition_2} - MA plot")
        y_max = max([abs(log2_fc) for log2_fc in deseq_result.log2FoldChange])
        plt.ylim(-1.1 * y_max, 1.1 * y_max)
        plt.xlabel("log10 base mean")
        plt.ylabel("log2 fold-change")
        self._pp_scatterplots.savefig()
        plt.close(fig)

    def create_volcano_plots(self, volcano_plot_path, volcano_plot_adj_path):
        self._pp_raw = PdfPages(volcano_plot_path)
        self._pp_adj = PdfPages(volcano_plot_adj_path)
        conditions = self._extract_condition_names()
        for condition_1 in conditions:
            for condition_2 in conditions:
                if condition_1 == condition_2:
                    continue
                self._create_volcano_plots(condition_1, condition_2)
        self._pp_raw.close()
        self._pp_adj.close()

    def _extract_condition_names(self):
        with open(self._deseq_script_path) as deseq_script_fh:
            for line in deseq_script_fh:
                if line.startswith("conds <- c('"):
                    for string in ["conds <- c('", "(", ")", "'", ","]:
                        line = line.replace(string, "")
                    return sorted(list(set(line.split())))

    def _create_volcano_plots(self, condition_1, condition_2):
        deseq_path = self._deseq_path_template + f"{condition_1}_vs_{condition_2}.csv"
        (
            basemean,
            log2_fold_changes,
            p_values,
            adj_p_values,
        ) = self._parse_deseq_file(deseq_path)
        cleaned_p_values = []
        cleaned_log2_fold_changes = []
        for p_value, log2_fold_change in zip(p_values, log2_fold_changes):
            if p_value != "NA" and log2_fold_change != "NA":
                cleaned_p_values.append(p_value)
                cleaned_log2_fold_changes.append(log2_fold_change)
        self._create_volcano_plot(
            cleaned_log2_fold_changes,
            cleaned_p_values,
            condition_1,
            condition_2,
            self._pp_raw,
        )
        cleaned_adj_p_values = []
        cleaned_log2_fold_changes = []
        for p_value, log2_fold_change in zip(adj_p_values, log2_fold_changes):
            if p_value != "NA" and log2_fold_change != "NA":
                cleaned_adj_p_values.append(p_value)
                cleaned_log2_fold_changes.append(log2_fold_change)
        self._create_volcano_plot(
            cleaned_log2_fold_changes,
            cleaned_adj_p_values,
            condition_1,
            condition_2,
            self._pp_adj,
            pvalue_string_mod="(adjusted)",
        )

    def _create_volcano_plot(
        self,
        log2_fold_changes,
        p_values,
        condition_1,
        condition_2,
        pp,
        pvalue_string_mod="",
    ):
        fig = plt.figure()
        max_log_2_fold_change = max(
            [abs(min(log2_fold_changes)), abs(max(log2_fold_changes))]
        )
        # To avoid problem with zero in log10
        p_values = np.array(p_values)
        p_values[p_values == 0.0] = 10**-100
        mod_p_values = -1 * np.log10(p_values)
        mod_p_values[mod_p_values == float("+inf")] = 0.0
        mod_p_values[mod_p_values == float("-inf")] = 0.0
        max_mod_p_values = max(mod_p_values)
        # Set axis ranges
        plt.axis(
            [
                -1 * max_log_2_fold_change,
                max_log_2_fold_change,
                0,
                max_mod_p_values,
            ]
        )
        plt.plot(log2_fold_changes, mod_p_values, "k.", alpha=0.3)
        # Add axis labels
        plt.xlabel("log$_2$ fold change")
        plt.ylabel(f"- log$_{10}$ p-value {pvalue_string_mod}")
        plt.title(f"{condition_1} vs. {condition_2}")
        significant_p_value = -1 * np.log10(self._p_value_significance_limit)
        plt.plot(
            [-1 * max_log_2_fold_change, max_log_2_fold_change],
            [significant_p_value, significant_p_value],
            linestyle="dotted",
            color="green",
            alpha=0.5,
        )
        # log2_fold_change_limit = 1 (never used)
        plt.plot(
            [
                -1 * self._log_2_fold_chance_limit,
                -1 * self._log_2_fold_chance_limit,
            ],
            [0, max_mod_p_values],
            linestyle="dotted",
            color="green",
            alpha=0.5,
        )
        plt.plot(
            [self._log_2_fold_chance_limit, self._log_2_fold_chance_limit],
            [0, max_mod_p_values],
            linestyle="dotted",
            color="green",
            alpha=0.5,
        )
        pp.savefig()
        plt.close(fig)

    def _parse_deseq_file(self, deseq_path):
        basemeans = []
        log_2_fold_changes = []
        p_values = []
        adj_p_values = []
        for row in csv.reader(open(deseq_path), delimiter="\t"):
            if row[0].startswith("baseMean"):
                continue
            basemeans.append(float(row[self._basemean_column]))
            for value_list, column_no in (
                (log_2_fold_changes, self._log_2_fold_chance_column),
                (p_values, self._p_value_column),
                (adj_p_values, self._adj_p_value_column),
            ):
                try:
                    value_list.append(float(row[column_no]))
                except ValueError:
                    value_list.append(row[column_no])
        return basemeans, log_2_fold_changes, p_values, adj_p_values
