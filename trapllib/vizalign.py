import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class AlignViz(object):

    def __init__(self, lib_names, read_processing_stats_path, 
                 read_aligner_stats_path):
        self._lib_names = lib_names
        self._read_processing_stats_path = read_processing_stats_path
        self._read_aligner_stats_path = read_aligner_stats_path

    def read_stat_files(self):
        with open(self._read_processing_stats_path) as read_processing_stats_fh:
            self._read_processing_stats = json.loads(
                read_processing_stats_fh.read())
        with open(self._read_aligner_stats_path) as read_aligner_stats_fh:
            self._read_aligner_stats = json.loads(
                read_aligner_stats_fh.read())

    def plot_input_read_length(self, plot_path):
        self._plot_read_lengths(
            "read_length_before_processing_and_freq", 
            "Length distribution of the input reads - %s",
            plot_path)

    def plot_processed_read_length(self, plot_path):
        self._plot_read_lengths(
            "read_length_after_processing_and_freq", 
            "Length distribution of the processed reads - %s",
            plot_path)

    def _plot_read_lengths(self, dict_key, title_template, output_file):
        pp = PdfPages(output_file)
        for lib in self._lib_names:
            lengths_and_freqs = self._read_processing_stats[lib][dict_key]
            fig = self._generate_histogram(lengths_and_freqs, title_template % lib)
            pp.savefig()
        pp.close()
    
    def _generate_histogram(self, lengths_and_freqs, title):
        fig = plt.figure()
        lengths = np.array([int(length) for length in lengths_and_freqs.keys()])
        freqs = np.array([int(freq) for freq in lengths_and_freqs.values()])
        ax = fig.add_subplot(111)
        plt.title(title)
        plt.xlabel("Read length [nt]")
        plt.ylabel("Frequency")
        font = {'size' : 8}
        matplotlib.rc('font', **font)
        ax.xaxis.set_ticks_position("bottom")
        plt.xticks(np.arange(0, max(lengths)+1, 10.0))
        plt.bar(lengths, freqs, align="center", color="black", edgecolor="none")
        plt.xlim([0,max(lengths)+1])
