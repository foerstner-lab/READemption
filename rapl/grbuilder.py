import concurrent.futures
from subprocess import call
from libs.sam import SamParser
from rapl.helper import Helper
from rapl.parameters import Parameters
from rapl.paths import Paths
from rapl.rapl_tools.sam2gr import Sam2Gr

class GrBuilder(object):

    def __init__(self):
        self.paths = Paths()
        self.parameters = Parameters()
        self.helper = Helper()
    
    def build_gr_files(self):
        """Generate GR files for all read/genome file combinations."""
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.parameters.python_number_of_threads) as executor:
            for read_file in self.paths.read_files:
                for genome_file in self.paths.genome_files:
                    executor.submit(self._build_gr_file, read_file, genome_file)

    def build_read_normalized_gr_files(self):
        """Generate normalized GR files for all read/genome files"""
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.parameters.python_number_of_threads) as executor:
            for genome_file in self.paths.genome_files:
                lowest_number_of_mappings = self.helper._lowest_number_of_mappings(
                    genome_file)
                for read_file in self.paths.read_files:
                    executor.submit(
                        self._build_read_normalized_gr_file, read_file, 
                        genome_file, lowest_number_of_mappings)

    def _build_gr_file(self, read_file, genome_file):
        """Generate GR files

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the first mapping file.
        - `genome_file`: name of the target genome file.
        """
        genome_file_header = self.helper.get_header_of_genome_file(genome_file)
        sam2gr = Sam2Gr(self.paths.combined_mapping_file_a_filtered(read_file),
                        mapping_target=genome_file_header,
                        output_prefix=self.paths.gr_file(read_file, genome_file))
        sam2gr.check_parameters()
        sam2gr.collect_intensities()
        sam2gr.print_output()
            
    def _build_read_normalized_gr_file(self, read_file, genome_file, 
                                  lowest_number_of_mappings):
        """Generate read normalized GR files

        Arguments:
        - `read_file`: orignal read file used to generate the mappings.
        - `genome_file,`: target genome file
        - `lowest_number_of_mappings`: the lowester number of mappings
                                       found for all read libs for a
                                       the genome file.

        """
        genome_file_header = self.helper.get_header_of_genome_file(genome_file)
        sam2gr = Sam2Gr(
            self.paths.combined_mapping_file_a_filtered(read_file),
            normalize_by_reads=True, 
            normalization_multiplier=lowest_number_of_mappings,
            mapping_target=genome_file_header, 
            output_prefix=self.paths.gr_read_normalized_file(
                read_file, genome_file))
        sam2gr.check_parameters()
        sam2gr.collect_intensities()
        sam2gr.print_output()
