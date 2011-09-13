import sys
sys.path.append(".")
from segemehl import Segemehl

class ReadMapper(object):

    def __init__(self):
        self.segemehl = Segemehl()
    
    def build_index(self, genome_file_paths, index_file_path):
        """
        
        Arguments:
        - `genome_file_paths`:
        - `index_file_path`:
        """
        self.segemehl.build_index(genome_file_paths, index_file_path)
        
    def run_mappings(self, read_file_paths, genome_file_paths, 
                     index_file_path, output_file_paths, nomatch_file_paths,
                     treads, accuracy):
        for read_file_path, output_file_path, nomatch_file_path in zip(
            read_file_paths, output_file_paths, nomatch_file_paths)
            self.segemehl.map_reads(
                read_file_path, index_file_path, genome_file_paths, 
                output_file_path, nonmatch_file=nomatch_file_path, 
                treads=treads, accuracy=accuracy)
