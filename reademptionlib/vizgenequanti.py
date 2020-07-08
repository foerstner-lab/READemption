import csv
from collections import defaultdict
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as colors
import matplotlib.cm as cm


class GeneQuantiViz(object):
    def __init__(
        self,
        gene_wise_quanti_combined_path,
        lib_names,
        use_antisene=True,
        axis_min=None,
        axis_max=None,
    ):
        self._gene_wise_quanti_combined_path = gene_wise_quanti_combined_path
        self._lib_names = lib_names
        self._use_antisene = use_antisene
        self._axis_min = axis_min
        self._axis_max = axis_max

    def parse_input_table(self):
        self._lib_names_and_countings = defaultdict(list)
        # Dict of dict of dict:
        # lib name -> relative direction (sense/anti-sense)
        # -> annotation type (CDS, rRNA ..)
        self._lib_names_and_class_quanti = defaultdict(
            lambda: defaultdict(lambda: defaultdict(float))
        )
        for row in csv.reader(
            open(self._gene_wise_quanti_combined_path), delimiter="\t"
        ):
            if row[0].startswith("Orientation"):
                continue
            for index, cell in enumerate(row[10:]):
                self._lib_names_and_countings[self._lib_names[index]].append(
                    float(cell)
                )
                self._lib_names_and_class_quanti[self._lib_names[index]][
                    row[0]
                ][row[3]] += float(cell)

    def plot_correlations(self, plot_path):
        self._prepare_document(plot_path)
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
        ax.set_yscale("log")
        ax.set_xscale("log")
        # Draw line
        plt.plot(
            [self._axis_min, self._axis_max],
            [self._axis_min, self._axis_max],
            linestyle="solid",
            color="green",
            alpha=0.4,
        )
        # Calculate the Pearson correlation coefficient
        corr_coeff = np.corrcoef(
            self._lib_names_and_countings[lib_1],
            self._lib_names_and_countings[lib_2],
        )[0][1]
        # Set axis ranges
        plt.axis(
            [self._axis_min, self._axis_max, self._axis_min, self._axis_max]
        )
        plt.title("%s vs. %s\n(r = %s)" % (lib_1, lib_2, corr_coeff))
        plt.plot(
            self._lib_names_and_countings[lib_1],
            self._lib_names_and_countings[lib_2],
            "k.",
            alpha=0.2,
        )
        plt.xlabel("Expression %s" % lib_1)
        plt.ylabel("Expression %s" % lib_2)
        self._pp.savefig()
        plt.close(fig)

    def _set_axis_max(self):
        self._axis_max = max(
            [
                max(counting)
                for counting in self._lib_names_and_countings.values()
            ]
        )

    def plot_annotation_class_quantification(self, plot_path):
        all_classes_sorted = set()
        no_of_libs = len(self._lib_names)
        for directions in self._lib_names_and_class_quanti.values():
            for classes_and_counting in directions.values():
                for anno_class in classes_and_counting.keys():
                    all_classes_sorted.add(anno_class)
        all_classes_sorted = sorted(list(all_classes_sorted))
        bottom = np.array([0] * no_of_libs)
        fig = plt.figure()
        ax = plt.subplot(111)
        font = {"family": "sans-serif", "weight": "normal", "size": 6}
        matplotlib.rc("font", **font)
        plt.title("Number of reads per RNA classes")
        color_map = plt.get_cmap("Set3")
        cNorm = colors.Normalize(
            vmin=0,
            vmax=(
                len(all_classes_sorted)
                * len(
                    self._lib_names_and_class_quanti[self._lib_names[0]].keys()
                )
            )
            - 1,
        )
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=color_map)
        color_index = 0
        for direction in self._lib_names_and_class_quanti[
            self._lib_names[0]
        ].keys():
            for anno_class in all_classes_sorted:
                countings = [
                    self._lib_names_and_class_quanti[lib][direction][anno_class]
                    for lib in self._lib_names
                ]
                color = scalarMap.to_rgba(color_index)
                plt.bar(
                    range(no_of_libs),
                    countings,
                    align="center",
                    bottom=bottom,
                    linewidth=0,
                    color=color,
                    width=0.5,
                    label=anno_class + " " + direction,
                )
                bottom = bottom + countings
                color_index += 1
        plt.xticks(
            np.array(range(no_of_libs)),
            self._lib_names,
            rotation=45,
            ha="right",
        )
        plt.tight_layout()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.xaxis.set_ticks_position("none")
        plt.legend(loc="upper right", frameon=False, ncol=4)
        fig.savefig(plot_path)
