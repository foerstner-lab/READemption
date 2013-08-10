import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class AlignViz(object):

    def __init__(self, lib_names, read_processing_stats_path, 
                 read_aligner_stats_path, input_read_length_plot_path,
                 processed_reads_length_plot_path):
        self._lib_names = lib_names
        self._read_processing_stats_path = read_processing_stats_path
        self._read_aligner_stats_path = read_aligner_stats_path
        self._input_read_length_plot_path = input_read_length_plot_path
        self._processed_read_length_plot_path = processed_reads_length_plot_path

    def read_stat_files(self):
        with open(self._read_processing_stats_path) as read_processing_stats_fh:
            self._read_processing_stats = json.loads(
                read_processing_stats_fh.read())
        with open(self._read_aligner_stats_path) as read_aligner_stats_fh:
            self._read_aligner_stats = json.loads(
                read_aligner_stats_fh.read())

    def plot_input_read_length(self):
        self._plot_read_lengths(
            "read_length_before_processing_and_freq", 
            "Length distribution of the input read - %s",
            self._input_read_length_plot_path)

    def plot_processed_read_length(self):
        self._plot_read_lengths(
            "read_length_after_processing_and_freq", 
            "Length distribution of the processed reads - %s",
            self._processed_read_length_plot_path)

    def _plot_read_lengths(self, dict_key, title_template, output_file):
        pp = PdfPages(output_file)
        for lib in self._lib_names:
            lengths_and_freqs = self._read_processing_stats[lib][dict_key]
            lengths =  []
            for length, freq in lengths_and_freqs.items():
                lengths = lengths + ([int(length)] * int(freq))
            fig = self._generate_histogram(lengths, title_template % lib)
            pp.savefig()
        pp.close()
    
    def _generate_histogram(self, list_of_values, title):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.title(title)
        font = {'size' : 6}
        matplotlib.rc('font', **font)
        hist, bins = np.histogram(
            list_of_values,bins=max(list_of_values)-1, 
            range=(0,max(list_of_values)-1))
        width = 0.7 * (bins[1]-bins[0])
        center = (bins[:-1] + bins[1:])/2
        plt.bar(center, hist, align = 'center', width = width, color="black")
