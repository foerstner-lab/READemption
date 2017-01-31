import csv
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as colors
import matplotlib.cm as cm
from bokeh.plotting import figure, ColumnDataSource
from bokeh.charts import Bar
from bokeh.layouts import column
from bokeh.models import (HoverTool, BoxZoomTool, ResetTool, PanTool,
                          WheelZoomTool)
from bokeh.resources import CDN
from bokeh.embed import file_html


class GeneQuantiViz(object):
    
    def __init__(self, gene_wise_quanti_combined_path, lib_names, output_path,
                 use_antisene=True, axis_min=None, axis_max=None):
        self._gene_wise_quanti_combined_path = gene_wise_quanti_combined_path
        self._output_path = output_path
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
            lambda: defaultdict(lambda: defaultdict(float)))
        for row in csv.reader(
                open(self._gene_wise_quanti_combined_path), delimiter="\t"):
            if row[0].startswith("Orientation"):
                continue
            for index, cell in enumerate(row[10:]):
                self._lib_names_and_countings[
                    self._lib_names[index]].append(float(cell))
                self._lib_names_and_class_quanti[
                    self._lib_names[index]][row[0]][row[3]] += float(cell)
        self.file_handle_bokeh()

    def file_handle_bokeh(self):
        data_bokeh_overview = []
        gene_quanti_combined_raw = pd.read_csv(
            self._gene_wise_quanti_combined_path, sep='\t', index_col=False)
        attr_dic = []
        for index, row in gene_quanti_combined_raw.iterrows():
            attr_dic.append(self._dictionary_attributes(row))
        df_attributes = pd.DataFrame(attr_dic)
        gene_quanti_combined_raw = pd.concat([
            gene_quanti_combined_raw, df_attributes], axis=1, join_axes=[
                gene_quanti_combined_raw.index])
        gene_quanti_combined_raw.drop("Attributes", axis=1, inplace=True)
        data_bokeh_overview.append(
            self._plotting_data_bokeh_overview(gene_quanti_combined_raw))
        self._plot_bokeh_read_no_per_lib(
            gene_quanti_combined_raw, data_bokeh_overview)

    def _plotting_data_bokeh_overview(self, gene_quanti_combined_raw):
        gene_quanti_combined_raw = gene_quanti_combined_raw.rename(
            columns={
                'Orientation of counted reads relative to the strand location '
                'of the annotation': 'Orientation'})
        orientation_list = []
        for orientation, df_rest in gene_quanti_combined_raw.groupby(
                'Orientation'):
            orientation_list.append(orientation)
        feature_list = []
        for feature, df_rest in gene_quanti_combined_raw.groupby('Feature'):
            feature_list.append(feature)
        library = []
        feature_orientation = []
        count = []
        for lib_name in self._lib_names:
            library.append(lib_name)
            for orientation in orientation_list:
                for feature in feature_list:
                    df_intermediate = gene_quanti_combined_raw[
                        (gene_quanti_combined_raw.Orientation == orientation) & (
                            gene_quanti_combined_raw.Feature == feature)]
                    feature_orientation.append('%s (%s)' % (
                        feature, orientation))
                    count.append(sum(df_intermediate[lib_name]))
        libraries = [
            item for item in library for i in range(len(feature_list) * len(
                orientation_list))]
        data_set = {}
        data_set['Library'] = libraries
        data_set['Features_Orientation'] = feature_orientation
        data_set['Gene count'] = count
        self.create_csv(data_set)
        self.plot_correlations()
        self.plot_annotation_class_quantification(orientation_list)
        return data_set

    def create_csv(self, data_set):
        with open(self._output_path + '/Read No per RNA class.csv',
                  'w') as csv_file:
            writer = csv.writer(csv_file, delimiter='\t')
            for key, value in data_set.items():
                writer.writerow([key] + value)
        
    def plot_correlations(self):
        self._prepare_document(self._output_path + '/Correlation.pdf')
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
        # Calculate the Pearson correlation coefficient
        corr_coeff = np.corrcoef(self._lib_names_and_countings[lib_1],
                                 self._lib_names_and_countings[lib_2])[0][1]
        # Set axis ranges
        plt.axis([self._axis_min, self._axis_max,
                  self._axis_min, self._axis_max])
        plt.title("%s vs. %s\n(r = %s)" % (lib_1, lib_2, corr_coeff))
        plt.plot(self._lib_names_and_countings[lib_1],
                 self._lib_names_and_countings[lib_2],
                 "k.", alpha=0.2)
        plt.xlabel("Expression %s" % lib_1)
        plt.ylabel("Expression %s" % lib_2)
        self._pp.savefig()
        plt.close(fig)

    def _set_axis_max(self):
        self._axis_max = max(
            [max(counting)
             for counting in self._lib_names_and_countings.values()])

    def plot_annotation_class_quantification(self, orientation_list):
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
        font = {'family': 'sans-serif', 'weight': 'normal', 'size': 6}
        matplotlib.rc('font', **font)
        plt.title("Number of reads per RNA classes")
        color_map = plt.get_cmap('Set3')
        cNorm = colors.Normalize(vmin=0, vmax=len(
            all_classes_sorted)*len(orientation_list)-1)
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=color_map)
        color_index = 0
        for direction in self._lib_names_and_class_quanti[
                self._lib_names[0]].keys():
            for anno_class in all_classes_sorted:
                countings = [
                    self._lib_names_and_class_quanti[lib][direction][
                        anno_class] for lib in self._lib_names]
                color = scalarMap.to_rgba(color_index)
                plt.bar(range(no_of_libs), countings, align="center",
                        bottom=bottom, linewidth=0, color=color, width=0.5,
                        label=anno_class+" "+direction)
                bottom = bottom + countings
                color_index += 1
        plt.xticks(np.array(range(no_of_libs)), self._lib_names, rotation=45,
                   ha="right")
        plt.tight_layout()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.xaxis.set_ticks_position("none")
        plt.legend(loc="upper right", frameon=False, ncol=4)
        fig.savefig(self._output_path + '/Read No per RNA Class.pdf')

    def _plot_bokeh_read_no_per_lib(
            self, gene_quanti_combined_raw, data_bokeh_overview):
        try:
            colors = ['Black', 'DarkSlateGray', 'MediumVioletRed', 'DarkCyan',
                      'Indigo', 'DarkBlue', '#a6cee3', '#1f78b4', '#b2df8a',
                      '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00',
                      '#cab2d6', '#6a3d9a', '#ffff99', '#b15928']
        except IndexError:
            sys.stderr.write("Appearently the number of features to be "
                             "plotted exceeds '18'. Exiting\n")
            sys.exit(2)
        plots = []
        for data_set in data_bokeh_overview:
            pl = Bar(data_set, values='Gene count', label='Library', agg='sum',
                     stack='Features_Orientation', palette=colors, tools=[
                         HoverTool(tooltips=[("No of reads", "@height")]),
                         BoxZoomTool(), ResetTool(), PanTool(),
                         WheelZoomTool()])
            pl.legend.background_fill_alpha = 0.5
            plots.append(pl)
        self._plot_bokeh_correlation(gene_quanti_combined_raw, plots)

    def _plot_bokeh_correlation(self, gene_quanti_combined_raw, plots):
        for lib_1 in self._lib_names:
            for lib_2 in self._lib_names:
                if lib_1 != lib_2:
                    corr_coeff = np.corrcoef(
                        gene_quanti_combined_raw[lib_1],
                        gene_quanti_combined_raw[lib_2])[0][1]
                    line_max = max(max(gene_quanti_combined_raw[lib_1]),
                                   max(gene_quanti_combined_raw[lib_2]))
                    pl = figure(title='%s vs %s (r=%s)' % (
                        lib_1, lib_2, corr_coeff), tools=[HoverTool(tooltips=[
                            ("Protein_ID", "@Name"),
                            ("Sequence type", "@gbkey"),
                            ("Product", "@product")]), PanTool(),
                                                          BoxZoomTool(),
                                                          WheelZoomTool(),
                                                          ResetTool()])
                    pl.scatter(gene_quanti_combined_raw[lib_1],
                               gene_quanti_combined_raw[lib_2], color='Black')
                    pl.xaxis.axis_label = 'Expression %s' % (lib_1)
                    pl.yaxis.axis_label = 'Expression %s' % (lib_2)
                    pl.title.text_font_style = 'italic'
                    pl.line([0.1, line_max], [0.1, line_max], line_width=0.2,
                            color='Green')
                    plots.append(pl)
        column(*plots)
        plot = column(*plots)
        html = file_html(plot, CDN, 'Viz Gene Quanti')
        with open(
                self._output_path + '/Visualization Gene Quantification.html',
                'w') as output_bokeh:
            output_bokeh.write(html)
    
    def _dictionary_attributes(self, row):
        dic = dict([key.split("=")
                    for key in row["Attributes"].split(";")])
        return dic

