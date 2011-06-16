import json
from subprocess import call
from rapl.grbuilder import GrBuilder
from rapl.parameters import Parameters
from rapl.paths import Paths
from rapl.helper import Helper
from libs.sam import SamParser
import concurrent.futures

class Annotations(object):

    def __init__(self):
        self.paths = Paths()
        self.parameters = Parameters()
        self.grbuilder = GrBuilder()
        self._get_annotation_files_from_config()

    def find_annotation_hits(self):
        """Search for overlaps of reads and annotations."""
        for read_file in self.paths.read_files:
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=self.parameters.python_number_of_threads) as executor:
                for annotation_file in self.annotation_files.keys():
                    executor.submit(
                        self._find_annotation_hits, read_file, annotation_file)

    def _find_annotation_hits(self, read_file, annotation_file):
        """Search for overlaps of reads and annotations.

        Arguments:
        - `read_file`: orignal read file used to generate the mappings.
        - `annotation_file`: an (NCBI) annotation file

        """
        genome_file_path = self.annotation_files[annotation_file]
        helper = Helper()
        genome_file_header = helper.get_header_of_genome_file(genome_file)
        call("%s %s/sam_hit_annotation_mapping.py -m %s -o %s %s %s \"%s\"" % (
                self.paths.python_bin, self.paths.bin_folder,
                self.parameters.min_overlap,
                self.paths.annotation_hit_file(read_file, annotation_file),
                self.paths.combined_mapping_file_a_filtered(read_file),
                self.paths.annotation_file(annotation_file),
                genome_file_header), shell=True)
    
    def _get_annotation_files_from_config(self):
        """Get the annations files from the config files.

        It extracts a dictionary that contains the names of the
        annotation files as keys and the names of the corresponding
        genome files as values.

        """
        self._read_config_file()
        self.annotation_files = self.config["annotation_and_genomes_files"]

    def _read_config_file(self):
        """Read the config file."""
        self.config = json.loads(open(self.paths.config_file).read())

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
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.parameters.python_number_of_threads) as executor:
            for annotation_file in self.annotation_files.keys():
                executor.submit(
                    self._build_annotation_hit_overview_read_normalized, 
                    annotation_file)

    def build_annotation_hit_overview_rpkm_normalized(self):
        """Create annotation hit overview tables RPKM normalized

        RPKM = Reads per kilobase of exon model per million mapped reads
        """
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.parameters.python_number_of_threads) as executor:
            for annotation_file in self.annotation_files.keys():
                executor.submit(
                    self._build_annotation_hit_overview_rpkm_normalized,
                    annotation_file)

    def _build_annotation_hit_overview(self, annotation_file):
        """Create annotation hit overview table. 

        Arguments:
        - `annotation_file`: an (NCBI) annotation file

        """
        annotation_hit_files_string = " ".join(
            [self.paths.annotation_hit_file(read_file, annotation_file) 
             for read_file in self.paths.read_files])
        call("%s %s/%s %s %s > %s" % (
                self.paths.python_bin, self.paths.bin_folder,
                "build_annotation_table_with_read_countings.py",
                self.paths.annotation_file(annotation_file),
                annotation_hit_files_string,
                self.paths.annotation_hit_overview_file(annotation_file)),
             shell=True)
        call("%s %s/%s -d a %s %s > %s" % (
                self.paths.python_bin, self.paths.bin_folder,
                "build_annotation_table_with_read_countings.py",
                self.paths.annotation_file(annotation_file),
                annotation_hit_files_string,
                self.paths.annotation_hit_overview_antisense_file(
                    annotation_file)), shell=True)

    def _build_annotation_hit_overview_read_normalized(self, annotation_file):
        """Create annotation hit overview table normalized by mapped
           reads.

        Arguments:
        - `annotation_file`: an (NCBI) annotation file

        """
        mapped_reads_counting_string = ":".join(
            [str(self.grbuilder._count_mapped_reads(read_file)) 
             for read_file in self.paths.read_files])
        annotation_hit_files_string = " ".join(
            [self.paths.annotation_hit_file(read_file, annotation_file) 
             for read_file in self.paths.read_files])
        call("%s %s/%s -n %s %s %s > %s" % (
                self.paths.python_bin, self.paths.bin_folder,
                "build_annotation_table_with_read_countings.py",
                mapped_reads_counting_string,
                self.paths.annotation_file(annotation_file),
                annotation_hit_files_string,
                self.paths.annotation_hit_overview_read_normalized_file(
                    annotation_file)), shell=True)
        call("%s %s/%s -d a -n %s %s %s > %s" % (
                self.paths.python_bin, self.paths.bin_folder,
                "build_annotation_table_with_read_countings.py",
                mapped_reads_counting_string,
                self.paths.annotation_file(annotation_file),
                annotation_hit_files_string,
                self.paths._annotation_hit_overview_read_normalized_antisense_file_path(
                    annotation_file)), shell=True)

    def _build_annotation_hit_overview_rpkm_normalized(self, annotation_file):
        """Create annotation hit overview table RPKM normalized. """
        mapped_reads_counting_string = ":".join(
            [str(self.grbuilder._count_mapped_reads(read_file)) 
             for read_file in self.paths.read_files])
        annotation_hit_files_string = " ".join(
            [self.paths.annotation_hit_file(read_file, annotation_file) 
             for read_file in self.paths.read_files])
        call("%s %s/%s -r -n %s %s %s > %s" % (
                self.paths.python_bin, self.paths.bin_folder,
                "build_annotation_table_with_read_countings.py",
                mapped_reads_counting_string,
                self.paths.annotation_file(annotation_file),
                annotation_hit_files_string,
                self.paths.annotation_hit_overview_rpkm_normalized_file(
                    annotation_file)),  shell=True)
        call("%s %s/%s -r -d a -n %s %s %s > %s" % (
                self.paths.python_bin, self.paths.bin_folder,
                "build_annotation_table_with_read_countings.py",
                mapped_reads_counting_string,
                self.paths.annotation_file(annotation_file),
                annotation_hit_files_string,
                self.paths.annotation_hit_overview_rpkm_normalized_antisense_file(
                    annotation_file)),   shell=True)

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
