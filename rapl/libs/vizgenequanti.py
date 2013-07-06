import numpy as np
from collections import defaultdict

class GeneQuantiViz(object):
    
    def __init__(self, gene_wise_quanti_combined_path, 
                 use_antisene=True, axis_min=None, axis_max=None):
        self._gene_wise_quanti_combined_path = gene_wise_quanti_combined_path
        self._use_antisene = use_antisene

    def parse_input_table(self):
        print(self._gene_wise_quanti_combined_path)
        self._lib_names_and_countings = defaultdict(np.array)
        self._feature_type = []

        
    def plot_correlations(self):
        pass

    def plot_feautere_distributions(self):
        pass
