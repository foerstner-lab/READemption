from subprocess import call
from libs.sam import SamParser
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
        call("%s %s/sam2gr.py -o %s %s" % (
                self.paths.python_bin, self.paths.bin_folder,
                self.paths.gr_file(read_file, genome_file),
                self.paths.combined_mapping_file_a_filtered_split(
                    read_file, genome_file)),
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
        call("%s %s/sam2gr.py -r -m %s -o %s %s" % (
                self.paths.python_bin, self.paths.bin_folder, lowest_number_of_mappings,
                self.paths.gr_read_normalized_file(read_file, genome_file),
                self.paths.combined_mapping_file_a_filtered_split(
                    read_file, genome_file)), shell=True)

    def build_nucl_normalized_gr_files(self):
        """Generate normalized GR files for all read/genome files"""
        for genome_file in self.paths.genome_files:
            lowest_number_of_mapped_nucleotides = (
                self._lowest_number_of_mapped_nucleotides(genome_file))
            for read_file in self.paths.read_files:
                self._build_nucl_normalized_gr_file(
                    read_file, genome_file, lowest_number_of_mapped_nucleotides)
            
    def _build_nucl_normalized_gr_file(self, read_file, genome_file, 
                                  lowest_number_of_mapped_nucleotides):
        """Generate read normalized GR files

        Arguments:
        - `read_file`: orignal read file used to generate the mappings.
        - `genome_file,`: target genome file
        - `lowest_number_of_mappings`: the lowester number of mappings
                                       found for all read libs for a
                                       the genome file.
        """
        call("%s %s/sam2gr.py -n -m %s -o %s %s" % (
                self.paths.python_bin, self.paths.bin_folder, 
                lowest_number_of_mapped_nucleotides,
                self.paths.gr_nucl_normalized_file(read_file, genome_file),
                self.paths.combined_mapping_file_a_filtered_split(
                    read_file, genome_file)), shell=True)

    def _lowest_number_of_mappings(self, genome_file):
        """Return the lowest number of mappings found.

        Arguments:
        - `genome_file`: target genome file

        """
        lowest_number_of_mappings = min(
            [self._count_mapped_reads(read_file, genome_file) 
             for read_file in self.paths.read_files])
        # Do avoid multiplication by zero
        if lowest_number_of_mappings == 0:
            lowest_number_of_mappings = 1
        return(lowest_number_of_mappings)

    def _lowest_number_of_mapped_nucleotides(self, genome_file):
        """Return the lowest number of mapping mapped nucleotides.

        Arguments:
        - `genome_file`: target genome file

        """
        lowest_number_of_mapped_nucleotides = min(
            [self._count_mapped_nucleotides(read_file, genome_file) 
             for read_file in self.paths.read_files])
        # Do avoid multiplication by zero
        if lowest_number_of_mapped_nucleotides == 0:
            lowest_number_of_nucleotides = 1
        return(lowest_number_of_mapped_nucleotides)

    def _count_mapped_reads(self, read_file, genome_file):
        """Count number of successfully mapped reads.

        Arguments:
        - `read_file`: orignal read file used to generate the mappings.
        - `genome_file`: targe genome file

        """
        sam_parser = SamParser()
        seen_ids = {}
        for entry in sam_parser.entries(
            self.paths.combined_mapping_file_a_filtered_split(
                read_file, genome_file)):
            seen_ids[entry['query']] = 1
        return(len(seen_ids))
    
    def _count_mapped_nucleotides(self, read_file, genome_file):
        """Count number of successfully mapped reads.

        Arguments:
        - `read_file`: orignal read file used to generate the mappings.
        - `genome_file`: targe genome file

        """
        sam_parser = SamParser()
        seen_ids = {}
        nucleotide_counting = 0
        for entry in sam_parser.entries(
            self.paths.combined_mapping_file_a_filtered_split(
                read_file, genome_file)):
            if entry['query'] in seen_ids:
                continue
            nucleotide_counting += len(entry['sequence'])
            seen_ids[entry['query']] = 1
        return(nucleotide_counting)

    def find_annotation_hits(self):
        """Search for overlaps of reads and annotations."""
        for read_file in self.paths.read_files:
            for annotation_file in self.annotation_files.keys():
                self._find_annotation_hits(read_file, annotation_file)
