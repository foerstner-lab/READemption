import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class DESeqViz(object):
    
    def __init__(self, deseq_script_path, deseq_path_template,
                 volcano_plot_file_path, use_antisene=True):
        self._deseq_script_path = deseq_script_path
        self._deseq_path_template = deseq_path_template
        self._volcano_plot_file_path = volcano_plot_file_path
        self._use_antisene = use_antisene
        self._log_2_fold_chance_column = 6
        self._p_value_column = 7
        self._p_adj_value_column = 8
        self._p_value_significance_limit = 0.05
        self._log_2_fold_chance_limit = 1

    def create_volcano_plots(self):
        self._pp = PdfPages(self._volcano_plot_file_path)
        conditions = self._extract_condition_names()
        for condition_1 in conditions:
            for condition_2 in conditions:
                if condition_1 == condition_2:
                    continue
                self._create_volcano_plot(condition_1, condition_2)
        self._pp.close()                

    def _extract_condition_names(self):
        with open(self._deseq_script_path) as deseq_script_fh:
            for line in deseq_script_fh:
                if line.startswith("conds <- c('"):
                    for string in ["conds <- c('", "(", ")", "'", ","]:
                        line = line.replace(string, "")
                    return sorted(list(set(line.split())))
        
    def _create_volcano_plot(self, condition_1, condition_2):
        deseq_path = self._deseq_path_template % (condition_1, condition_2)
        log2_fold_changes, p_values, p_adj_values = self._parse_deseq_file(
            deseq_path)
        max_log_2_fold_change = max(
            [abs(min(log2_fold_changes)), 
             abs(max(log2_fold_changes))])
        mod_p_values = -1 * np.log10(p_values)
        max_mod_p_values = max(mod_p_values)
        # Set axis ranges
        plt.axis([-1*max_log_2_fold_change, max_log_2_fold_change, 
                   0, max_mod_p_values])
        plt.plot(log2_fold_changes, mod_p_values, "k.", alpha=0.3)
        # Add axis labels
        plt.xlabel("log$_2$ fold change")
        plt.ylabel("- log$_{10}$ p-value")
        plt.title("%s vs. %s" % (condition_1, condition_2))
        signifant_p_value = -1*np.log10(self._p_value_significance_limit)
        plt.plot([-1*max_log_2_fold_change, max_log_2_fold_change], 
                 [signifant_p_value, signifant_p_value], 
                 linestyle="dotted", color="green", alpha=0.5)
        log2_fold_change_limit = 1
        plt.plot([-1*self._log_2_fold_chance_limit, 
                   -1*self._log_2_fold_chance_limit],
                 [0, max_mod_p_values], 
                 linestyle="dotted", color="green", alpha=0.5)
        plt.plot([self._log_2_fold_chance_limit, 
                  self._log_2_fold_chance_limit], 
                 [0, max_mod_p_values], 
                 linestyle="dotted", color="green", alpha=0.5)
        self._pp.savefig()

    def _parse_deseq_file(self, deseq_path):
        log_2_fold_changes = []
        p_value = []
        p_adj_value = []
        for row in csv.reader(open(deseq_path), delimiter="\t"):
            if row[0].startswith("id"):
                continue
            skip = False
            for column_no in [
                self._log_2_fold_chance_column, self._p_value_column,
                self._p_adj_value_column]:
                if row[column_no] in ["NA", "-Inf", "Inf"]:
                    skip = True
                    break
            if skip is True:
                continue
            log_2_fold_changes.append(float(row[self._log_2_fold_chance_column]))
            p_value.append(float(row[self._p_value_column]))
            p_adj_value.append(float(row[self._p_adj_value_column]))
        return log_2_fold_changes, p_value,  p_adj_value
