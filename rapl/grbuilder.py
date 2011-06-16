from subprocess import call
from libs.sam import SamParser
from rapl.helper import Helper
from rapl.parameters import Parameters
from rapl.paths import Paths

class GrBuilder(object):

    def __init__(self):
        self.paths = Paths()
    
    def build_gr_files(self):
        """Generate GR files for all read/genome file combinations."""
        for read_file in self.paths.read_files:
            for genome_file in self.paths.genome_files:
                self._build_gr_file(read_file, genome_file)

    def _build_gr_file(self, read_file, genome_file):
        """Generate GR files

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the first mapping file.
        - `genome_file`: name of the target genome file.
        """
        helper = Helper()
        genome_file_header = helper.get_header_of_genome_file(genome_file)
        call("%s %s/sam2gr.py -t %s -o %s %s" % (
                self.paths.python_bin, self.paths.bin_folder, 
                genome_file_header, self.paths.gr_file(read_file, genome_file),
                self.paths.combined_mapping_file_a_filtered(read_file)), 
             shell=True)

    def build_read_normalized_gr_files(self):
        """Generate normalized GR files for all read/genome files"""
        for genome_file in self.paths.genome_files:
            lowest_number_of_mappings = self._lowest_number_of_mappings(
                genome_file)
            for read_file in self.paths.read_files:
                self._build_read_normalized_gr_file(
                    read_file, genome_file, lowest_number_of_mappings)
            
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
        helper = Helper()
        genome_file_header = helper.get_header_of_genome_file(genome_file)
        call("%s %s/sam2gr.py -t %s -r -m %s -o %s %s" % (
                self.paths.python_bin, self.paths.bin_folder, genome_file_header,
                lowest_number_of_mappings, 
                self.paths.gr_read_normalized_file(read_file, genome_file),
                self.paths.combined_mapping_file_a_filtered(read_file)), shell=True)

    def _lowest_number_of_mappings(self, genome_file):
        """Return the lowest number of mappings found.

        Arguments:
        - `genome_file`: target genome file

        """
        lowest_number_of_mappings = min(
            [self._count_mapped_reads(read_file) 
             for read_file in self.paths.read_files])
        # Do avoid multiplication by zero
        if lowest_number_of_mappings == 0:
            lowest_number_of_mappings = 1
        return(lowest_number_of_mappings)

    def _count_mapped_reads(self, read_file):
        """Count number of successfully mapped reads.

        Arguments:
        - `read_file`: orignal read file used to generate the mappings.

        """
        sam_parser = SamParser()
        return(
            sam_parser.number_of_mapped_reads(
                self.paths.combined_mapping_file_a_filtered(read_file)))
    
    # def _count_mapped_nucleotides(self, read_file, genome_file):
    #     """Count number of successfully mapped reads.

    #     Arguments:
    #     - `read_file`: orignal read file used to generate the mappings.
    #     - `genome_file`: targe genome file

    #     """
    #     sam_parser = SamParser()
    #     nucleotide_counting = 0
    #     prev_entry = ""
    #     for entry in sam_parser.entries(
    #         self.paths.combined_mapping_file_a_filtered_split(
    #             read_file, genome_file)):
    #         if entry['query'] == prev_entry:
    #             continue
    #         nucleotide_counting += len(entry['sequence'])
    #         prev_entry = entry['query']
    #     return(nucleotide_counting)

    # def _lowest_number_of_mapped_nucleotides(self, genome_file):
    #     """Return the lowest number of mapping mapped nucleotides.

    #     Arguments:
    #     - `genome_file`: target genome file

    #     """
    #     lowest_number_of_mapped_nucleotides = min(
    #         [self._count_mapped_nucleotides(read_file, genome_file) 
    #          for read_file in self.paths.read_files])
    #     # Do avoid multiplication by zero
    #     if lowest_number_of_mapped_nucleotides == 0:
    #         lowest_number_of_nucleotides = 1
    #     return(lowest_number_of_mapped_nucleotides)

