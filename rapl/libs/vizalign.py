import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class AlignViz(object):

    def __init__(self, read_processing_stats_path, read_aligner_stats_path):
        self._read_processing_stats_path = read_processing_stats_path
        self._read_aligner_stats_path = read_aligner_stats_path

    def read_stat_files(self):
        self._read_processing_stats = None
        self._read_aligner_stats = None

    def plot_input_read_length(self):
        pass

    def plot_processed_read_length(self):
        pass
    
