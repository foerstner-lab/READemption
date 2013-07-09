import csv
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class GeneQuantiViz(object):
    
    def __init__(self, gene_wise_quanti_combined_path, lib_names, 
                 scatter_plot_file_path, use_antisene=True, 
                 axis_min=None, axis_max=None):
        self._gene_wise_quanti_combined_path = gene_wise_quanti_combined_path
        self._lib_names = lib_names
        self._scatter_plot_file_path = scatter_plot_file_path
        self._use_antisene = use_antisene
        self._axis_min = axis_min
        self._axis_max = axis_max

    def parse_input_table(self):
        self._lib_names_and_countings = defaultdict(list)
        for row in csv.reader(
            open(self._gene_wise_quanti_combined_path), delimiter="\t"):
            if len(row[0]) == 0:
                continue
            for index, cell in enumerate(row[10:]):
                self._lib_names_and_countings[
                    self._lib_names[index]].append(float(cell))
        
    def plot_correlations(self):
        self._prepare_document(self._scatter_plot_file_path)
        if self._axis_min is None:
            self._axis_min = 0.1
        if self._axis_max is None:
            self._set_axis_max()
        for lib_1 in self._lib_names:
            for lib_2 in self._lib_names:
                if lib_1 == lib_2:
                    continue
                self._plot_correlation(lib_1, lib_2)
        self._close_document()

    def _prepare_document(self, file_name):
        self._pp = PdfPages(file_name)

    def _close_document(self):
        self._pp.close()

    def _plot_correlation(self, lib_1, lib_2):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_aspect(1)
        ax.set_yscale('log')
        ax.set_xscale('log')
        # Draw line
        plt.plot([self._axis_min, self._axis_max], 
                 [self._axis_min, self._axis_max], 
                 linestyle="solid", color="green", alpha=0.4)
        # Set axis ranges
        plt.axis([self._axis_min, self._axis_max, 
                  self._axis_min, self._axis_max])
        plt.title("%s vs. %s" % (lib_1, lib_2))
        plt.plot(self._lib_names_and_countings[lib_1], 
                 self._lib_names_and_countings[lib_2], 
                 "k.", alpha=0.2)
        plt.xlabel("Expression %s" % lib_1)
        plt.ylabel("Expression %s" % lib_2)
        self._pp.savefig()

    def _set_axis_max(self):
        self._axis_max = max(
            [max(counting) 
             for counting in self._lib_names_and_countings.values()])
