import pandas as pd
import numpy as np
from bokeh.layouts import column
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import (HoverTool, BoxSelectTool, BoxZoomTool, ResetTool,
                          PanTool, WheelZoomTool, ResizeTool)
from bokeh.resources import CDN
from bokeh.embed import file_html
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


class DESeqViz(object):

    def __init__(self, input_file, output_path, cutoff, condition, alpha,
                 color_sig, color_non_sig, shape, glyph_size):
        self._input_file = input_file
        self._output_path = output_path
        self._cutoff = cutoff
        self._condition = condition
        self._alpha = alpha
        self._col_sig = str(color_sig)
        self._col_non_sig = str(color_non_sig)
        self._shape = shape
        self._glyph_size = glyph_size
        
    def read_and_modificate_input(self):
        deseq_data = pd.read_csv(self._input_file, sep='\t', skiprows=2)
        attr_dic = []
        for index, row in deseq_data.iterrows():
            attr_dic.append(self._dictionary_attributes(row))
        df_attributes = pd.DataFrame(attr_dic)
        deseq_data_mod = pd.concat([deseq_data, df_attributes],
                                   axis=1, join_axes=[deseq_data.index])
        deseq_data_mod.drop("Attributes", axis=1, inplace=True)
        self._plot_MH_matplotlib(deseq_data_mod)
        self._plot_MA_matplotlib(deseq_data_mod)
        self._plot_bokeh_MA(deseq_data_mod)
    
    def _plot_bokeh_MA(self, deseq_data_mod):
        plots = []
        deseq_sig = deseq_data_mod[deseq_data_mod.padj < self._cutoff]
        deseq_no_sig = deseq_data_mod[deseq_data_mod.padj >= self._cutoff]
        pl = figure(tools=[HoverTool(tooltips=[
            ("Protein_ID", "@Name"), ("Sequence type", "@gbkey"),
            ("Product", "@product"), ("log2FoldChange", "@log2FoldChange"),
            ("baseMean", "@baseMean")]), PanTool(), BoxSelectTool(),
                           BoxZoomTool(), WheelZoomTool(), ResizeTool(),
                           ResetTool()])
        pl.background_fill_color = "White"
        pl.grid.grid_line_color = "black"
        pl.xaxis.axis_label = 'log10 baseMean'
        pl.yaxis.axis_label = 'log2 FoldChange'
        pl.title.text = 'MA Plot'
        self._plot(pl, np.log10(deseq_sig["baseMean"]),
                   deseq_sig["log2FoldChange"], self._col_sig,
                   self._glyph_size, self._alpha, 'padj significant',
                   ColumnDataSource(deseq_sig))
        self._plot(pl, np.log10(deseq_no_sig["baseMean"]),
                   deseq_no_sig["log2FoldChange"], self._col_non_sig,
                   self._glyph_size, self._alpha, 'padj non-significant',
                   ColumnDataSource(deseq_no_sig))
        plots.append(pl)
        self._plot_bokeh_MH(deseq_data_mod, plots)

    def _plot_bokeh_MH(self, deseq_data_mod, plots):
        for replicon, data_group in deseq_data_mod.groupby(["Sequence name"]):
            data_group_sig = data_group[data_group.padj < self._cutoff]
            data_group_no_sig = data_group[data_group.padj >= self._cutoff]
            pl = figure(tools=[HoverTool(tooltips=[
                ("Protein_ID", "@Name"), ("Sequence type", "@gbkey"),
                ("Product", "@product"), ("log2FoldChange", "@log2FoldChange"),
                ("baseMean", "@baseMean")]), PanTool(), BoxSelectTool(),
                               BoxZoomTool(), WheelZoomTool(), ResizeTool(),
                               ResetTool()])
            pl.background_fill_color = "White"
            pl.grid.grid_line_color = "black"
            pl.xaxis.axis_label = 'Sequence Start Position'
            pl.yaxis.axis_label = 'log2 FoldChange'
            pl.title.text = 'MH Plot (' + replicon + ')'
            pl.circle(data_group_sig["Start"],
                      data_group_sig["log2FoldChange"], alpha=float(0.5),
                      size=self._calc_glyph_size(data_group_sig),
                      legend=('padj significant (cutoff: ' + str(
                          self._cutoff) + ')'), color='Red',
                      source=ColumnDataSource(data_group_sig))
            pl.circle(data_group_no_sig["Start"],
                      data_group_no_sig["log2FoldChange"],
                      size=self._calc_glyph_size(
                          data_group_no_sig), alpha=float(0.5),
                      legend=('padj non-significant (cutoff: ' + str(
                          self._cutoff) + ')'), color='Black',
                      source=ColumnDataSource(data_group_no_sig))
            plots.append(pl)
        plot = column(*plots)
        html = file_html(plot, CDN, 'MA & MH Plot {}'.format(
            self._condition))
        with open('{}/MA & MH Plot {}.html'.format(self._output_path,
                                                   self._condition),
                  'w') as output_bokeh:
            output_bokeh.write(html)
        
    def _plot(self, pl, x, y, color, size, alpha, legend, source):
        if self._shape == 'circle':
            return(pl.circle(x, y, color=color, size=size, alpha=alpha,
                             legend=legend, source=source))
        if self._shape == 'square':
            return(pl.square(x, y, color=color, size=size, alpha=alpha,
                             legend=legend, source=source))

    def _dictionary_attributes(self, row):
        dic = dict([key.split("=")
                    for key in row["Attributes"].split(";")])
        return dic

    def _calc_glyph_size(self, data):
        return np.log2(data["End"].sub(data["Start"], axis=0))

    def _plot_MH_matplotlib(self, deseq_data_mod):
        with PdfPages('{}/MH_plot {}.pdf'.format(self._output_path,
                                                 self._condition)) as pdf:
            for replicon, data_group in deseq_data_mod.groupby([
                    "Sequence name"]):
                data_group_sig = data_group[data_group.padj < self._cutoff]
                data_group_no_sig = data_group[data_group.padj >= self._cutoff]
                matplotlib.style.use('ggplot')
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.scatter(data_group_sig["Start"].tolist(),
                           data_group_sig["log2FoldChange"].tolist(),
                           s=8, c='r', label=(
                               'padj significant(cutoff: ' + str(
                                   self._cutoff) + ')'))
                ax.scatter(data_group_no_sig["Start"].tolist(),
                           data_group_no_sig["log2FoldChange"].tolist(),
                           s=8, c='k', label=(
                               'padj non-significant(cutoff: ' + str(
                                   self._cutoff) + ')'))
                plt.xlabel('Sequence Start Position')
                plt.ylabel('log2FoldChange')
                plt.title('MH plot (' + replicon + ')')
                leg = plt.legend(fancybox=True, loc='best')
                leg.get_frame().set_alpha(0.5)
                pdf.savefig()
                plt.close()

    def _plot_MA_matplotlib(self, deseq_data_mod):
        deseq_sig = deseq_data_mod[deseq_data_mod.padj < self._cutoff]
        deseq_no_sig = deseq_data_mod[deseq_data_mod.padj >= self._cutoff]
        with PdfPages('{}/MA_plot {}.pdf'.format(self._output_path,
                                                 self._condition)) as pdf:
            matplotlib.style.use('ggplot')
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.scatter(np.log10(deseq_sig["baseMean"].tolist()),
                       deseq_sig["log2FoldChange"].tolist(), s=8, c='r',
                       label=('padj significant(cutoff: ' + str(
                           self._cutoff) + ')'))
            ax.scatter(np.log10(deseq_no_sig["baseMean"].tolist()),
                       deseq_no_sig["log2FoldChange"], s=8, c='k', label=(
                           'padj non-significant(cutoff: ' + str(
                               self._cutoff) + ')'))
            plt.xlabel('log10 baseMean')
            plt.ylabel('log2 FoldChange')
            plt.title('MA plot')
            leg = plt.legend(fancybox=True, loc='best')
            leg.get_frame().set_alpha(0.5)
            pdf.savefig()
            plt.close()
