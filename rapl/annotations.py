import concurrent.futures
import json
import sys
from subprocess import call
from rapl.grbuilder import GrBuilder
from rapl.parameters import Parameters
from rapl.paths import Paths
from rapl.helper import Helper
from libs.sam import SamParser
from rapl.rapl_tools.build_annotation_table_with_read_countings import AnnotationMappingTableBuilder

class Annotations(object):

    def __init__(self):
        self.paths = Paths()
        self.parameters = Parameters()
        self.grbuilder = GrBuilder()
        self.helper = Helper()
        self._get_annotation_files_from_config()

    def find_annotation_hits(self):
        """Search for overlaps of reads and annotations."""
        # Create a thread for each read file / annotation file combo.
        threads = []
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.parameters.python_number_of_threads) as executor:
            for read_file in self.paths.read_files:
                for annotation_file in self.annotation_files.keys():
                    threads.append(executor.submit(
                        self._find_annotation_hits, read_file, annotation_file))
        # Evaluate thread outcome
        self.helper.check_thread_completeness(threads)

    def build_annotation_hit_overview(self):
        """Create annotation hit overview tables."""
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.parameters.python_number_of_threads) as executor:
            for annotation_file in self.annotation_files.keys():
                executor.submit(
                    self._build_annotation_hit_overview, annotation_file)

    def build_annotation_hit_overview_read_normalized(self):
        """Create annotation hit overview tables normalized by mapped
           reads.

        """
        # Create a thread for each read file
        threads = []
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.parameters.python_number_of_threads) as executor:
            for annotation_file in self.annotation_files.keys():
                threads.append(executor.submit(
                    self._build_annotation_hit_overview_read_normalized, 
                    annotation_file))
        # Evaluate thread outcome
        self.helper.check_thread_completeness(threads)

    def build_annotation_hit_overview_nucleotide_normalized(self):
        """Create annotation hit overview tables - nucleotide based
        
        Overlapping nucleotides are counted and normalized by the
        total number of mapped nucleotides.
        """
        # Create a thread for each read file
        threads = []
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.parameters.python_number_of_threads) as executor:
            for annotation_file in self.annotation_files.keys():
                threads.append(executor.submit(
                    self._build_annotation_hit_overview_nucleotide_normalized, 
                    annotation_file))
        # Evalutate thread outcome
        self.helper.check_thread_completeness(threads)

    def build_annotation_hit_overview_rpkm_normalized(self):
        """Create annotation hit overview tables RPKM normalized

        RPKM = Reads per kilobase of exon model per million mapped reads
        """
        # Create a thread for each read file
        threads = []
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.parameters.python_number_of_threads) as executor:
            for annotation_file in self.annotation_files.keys():
                threads.append(executor.submit(
                    self._build_annotation_hit_overview_rpkm_normalized,
                    annotation_file))
        # Evalutate thread outcome
        self.helper.check_thread_completeness(threads)

    def _find_annotation_hits(self, read_file, annotation_file):
        """Search for overlaps of reads and annotations.

        Arguments:
        - `read_file`: orignal read file used to generate the mappings.
        - `annotation_file`: an (NCBI) annotation file

        """
        genome_file = self.annotation_files[annotation_file]
        genome_file_header = self.helper.get_header_of_genome_file(genome_file)
        call("%s %s/sam_hit_annotation_mapping.py -m %s -o %s %s %s \"%s\"" % (
                self.paths.python_bin, self.paths.bin_folder,
                self.parameters.min_overlap,
                self.paths.annotation_hit_file(read_file, annotation_file),
                self.paths.combined_mapping_file_a_filtered(read_file),
                self.paths.annotation_file(annotation_file),
                genome_file_header), shell=True)
    
    def _get_annotation_files_from_config(self):
        """Get the annotations files from the config files.

        It extracts a dictionary that contains the names of the
        annotation files as keys and the names of the corresponding
        genome files as values.

        """
        self._read_config_file()
        self.annotation_files = self.config["annotation_and_genomes_files"]

    def _read_config_file(self):
        """Read the config file."""
        self.config = json.loads(open(self.paths.config_file).read())

    def _build_annotation_hit_overview(self, annotation_file):
        """Create annotation hit overview table. 

        Arguments:
        - `annotation_file`: an (NCBI) annotation file

        """
        annotation_hit_files = [
            self.paths.annotation_hit_file(read_file, annotation_file) 
            for read_file in self.paths.read_files]
        # Sense
        annotation_table_builder = AnnotationMappingTableBuilder(
            self.paths.annotation_file(annotation_file), annotation_hit_files,
             strand_orientation="s",
            output_file=self.paths.annotation_hit_overview_file(annotation_file))
        annotation_table_builder.read_annotation_mapping_files()
        annotation_table_builder.read_annotation_file_and_print_output()
        # Antisense
        annotation_table_builder = AnnotationMappingTableBuilder(
            self.paths.annotation_file(annotation_file), annotation_hit_files,
            strand_orientation="a", 
            output_file=self.paths.annotation_hit_overview_antisense_file(
                annotation_file))
        annotation_table_builder.read_annotation_mapping_files()
        annotation_table_builder.read_annotation_file_and_print_output()

    def _build_annotation_hit_overview_read_normalized(self, annotation_file):
        """Create annotation hit overview table normalized by mapped
           reads.
n
        Arguments:
        - `annotation_file`: an (NCBI) annotation file

        """

        annotation_hit_files = [
            self.paths.annotation_hit_file(read_file, annotation_file) 
            for read_file in self.paths.read_files]
        mapped_reads_countings = [
            self.grbuilder._count_mapped_reads(read_file) 
            for read_file in self.paths.read_files]
        # Sense
        annotation_table_builder = AnnotationMappingTableBuilder(
            self.paths.annotation_file(annotation_file), 
            annotation_hit_files, strand_orientation="s", 
            normalization_factors=mapped_reads_countings,
            output_file=self.paths.annotation_hit_overview_read_normalized_file(
                annotation_file))
        annotation_table_builder.read_annotation_mapping_files()
        annotation_table_builder.read_annotation_file_and_print_output()
        # Antisense
        annotation_table_builder = AnnotationMappingTableBuilder(
            self.paths.annotation_file(annotation_file), annotation_hit_files,
            strand_orientation="a", normalization_factors=mapped_reads_countings,
            output_file=self.paths._annotation_hit_overview_read_normalized_antisense_file_path(
                     annotation_file))
        annotation_table_builder.read_annotation_mapping_files()
        annotation_table_builder.read_annotation_file_and_print_output()

    def _build_annotation_hit_overview_nucleotide_normalized(self, annotation_file):
        """Create annotation hit overview table - nucleotide based

        Arguments:
        - `annotation_file`: an (NCBI) annotation file

        """
        annotation_hit_files = [
            self.paths.annotation_hit_file(read_file, annotation_file) 
            for read_file in self.paths.read_files]
        mapped_nucleotide_countings = [
            self.helper.count_mapped_nucleotides(read_file) 
            for read_file in self.paths.read_files]
        # Sense
        annotation_table_builder = AnnotationMappingTableBuilder(
            self.paths.annotation_file(annotation_file), 
            annotation_hit_files, strand_orientation="s", 
            count_nucleotides=True,
            normalization_factors=mapped_nucleotide_countings,
            output_file=self.paths.annotation_hit_overview_nucl_normalized_file(
                annotation_file))
        annotation_table_builder.read_annotation_mapping_files()
        annotation_table_builder.read_annotation_file_and_print_output()
        # Antisense
        annotation_table_builder = AnnotationMappingTableBuilder(
            self.paths.annotation_file(annotation_file), annotation_hit_files,
            strand_orientation="a", normalization_factors=mapped_nucleotide_countings,
            count_nucleotides=True,
            output_file=self.paths.annotation_hit_overview_nucl_normalized_antisense_file(
                     annotation_file))
        annotation_table_builder.read_annotation_mapping_files()
        annotation_table_builder.read_annotation_file_and_print_output()

    def _build_annotation_hit_overview_rpkm_normalized(self, annotation_file):
        """Create annotation hit overview table RPKM normalized. """
        annotation_hit_files = [
            self.paths.annotation_hit_file(read_file, annotation_file) 
            for read_file in self.paths.read_files]
        mapped_reads_countings = [
            self.grbuilder._count_mapped_reads(read_file) 
            for read_file in self.paths.read_files]
        # Sense
        annotation_table_builder = AnnotationMappingTableBuilder(
            self.paths.annotation_file(annotation_file), 
            annotation_hit_files, strand_orientation="s", rpkm=True,
            normalization_factors=mapped_reads_countings, 
            output_file=self.paths.annotation_hit_overview_rpkm_normalized_file(
                annotation_file))
        annotation_table_builder.read_annotation_mapping_files()
        annotation_table_builder.read_annotation_file_and_print_output()
        # Antisense
        annotation_table_builder = AnnotationMappingTableBuilder(
            self.paths.annotation_file(annotation_file), annotation_hit_files,
            strand_orientation="a", rpkm=True,
            normalization_factors=mapped_reads_countings,
            output_file=self.paths.annotation_hit_overview_rpkm_normalized_antisense_file(
                annotation_file))
        annotation_table_builder.read_annotation_mapping_files()
        annotation_table_builder.read_annotation_file_and_print_output()

    def count_reads_in_intergenic_regions(self):
        """Get the number of reads in intergenic region (IGR)

        These are the reads that do not overlap with any annotations.
        """
        for read_file in self.paths.read_files:
            self._count_reads_in_intergenic_region(read_file)

    def _count_reads_in_intergenic_region(self, read_file):
        # TODO
        genome_and_annotation_files = self._group_annotation_files_by_genome_file()
        for genome_file, annotation_files in genome_and_annotation_files.items():
            mapped_reads = self._get_mapped_reads(read_file, genome_file)
            reads_with_annoation_overlaps = self._get_reads_with_annotation_overlap(
                read_file, annotation_files)
            reads_without_annotation_overlap = filter(
                lambda read: read in reads_with_annoation_overlaps, 
                mapped_reads.keys())
            fh = open(self.paths.no_annotation_hit_file(read_file, genome_file), "w")
            fh.write("\n".join(reads_without_annotation_overlap) + "\n")

    def _group_annotation_files_by_genome_file(self):
        genome_and_annotation_files = {}
        for annotation_file, genom_file in self.annotation_files.items():
            genome_and_annotation_files.setdefault(genom_file, [])
            genome_and_annotation_files[genom_file].append(annotation_file)
        return(genome_and_annotation_files)

    def _get_mapped_reads(self, read_file, genome_file):
        sam_parser = SamParser()
        mapped_reads = {}
        for entry in sam_parser.entries(
            self.paths.combined_mapping_file_a_filtered_split(
                read_file, genome_file)):
            mapped_reads[entry["query"]] = 1
        return(mapped_reads)
            
    def _get_reads_with_annotation_overlap(
        self, read_file, annotation_files):
        reads_with_annotation = {}
        for annotation_file in annotation_files:
            for line in open(self.paths.annotation_hit_file(
                    read_file, annotation_file)):
                if line[0] in ["#", "\n"]: continue
                reads_with_annotation[line.split("\t")[0]] = 1
        return(reads_with_annotation)
