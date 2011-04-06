import json
from subprocess import call
from rapl.grbuilder import GrBuilder
from rapl.parameters import Parameters
from rapl.pathes import Pathes

class Annotations(object):

    def __init__(self):
        self.pathes = Pathes()
        self.parameters = Parameters()
        self.grbuilder = GrBuilder()
        self._get_annotation_files_from_config()

    def find_annotation_hits(self):
        """Search for overlaps of reads and annotations."""
        for read_file in self.pathes.read_files:
            for annotation_file in self.annotation_files.keys():
                self._find_annotation_hits(read_file, annotation_file)

    def _find_annotation_hits(self, read_file, annotation_file):
        """Search for overlaps of reads and annotations.

        Arguments:
        - `read_file`: orignal read file used to generate the mappings.
        - `annotation_file`: an (NCBI) annotation file

        """
        genome_file = self.annotation_files[annotation_file]
        call("%s %s/segemehl_hit_annotation_mapping.py -m %s -o %s %s %s" % (
                self.pathes.python_bin, self.pathes.bin_folder,
                self.parameters.min_overlap,
                self.pathes.annotation_hit_file(read_file, annotation_file),
                self.pathes.combined_mapping_file_a_filtered_split(
                    read_file, genome_file),
                self.pathes.annotation_file(annotation_file)), shell=True)

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
        self.config = json.loads(open(self.pathes.config_file).read())

    def build_annotation_hit_overview(self):
        """Create annotation hit overview tables."""
        for annotation_file in self.annotation_files.keys():
            self._build_annotation_hit_overview(annotation_file)

    def build_annotation_hit_overview_read_normalized(self):
        """Create annotation hit overview tables normalized by mapped
           reads.

        """
        for annotation_file in self.annotation_files.keys():
            self._build_annotation_hit_overview_read_normalized(annotation_file)

    def build_annotation_hit_overview_nucl_normalized(self):
        """Create annotation hit overview tables normalized by mapped
        nucleotides.

        """
        for annotation_file in self.annotation_files.keys():
            self._build_annotation_hit_overview_nucl_normalized(annotation_file)

    def _build_annotation_hit_overview(self, annotation_file):
        """Create annotation hit overview table. 

        Arguments:
        - `annotation_file`: an (NCBI) annotation file

        """
        genome_file = self.annotation_files[annotation_file]
        annotation_hit_files_string = " ".join(
            [self.pathes.annotation_hit_file(read_file, annotation_file) 
             for read_file in self.pathes.read_files])
        call("%s %s/%s %s %s > %s" % (
                self.pathes.python_bin, self.pathes.bin_folder,
                "build_annotation_table_with_read_countings.py",
                self.pathes.annotation_file(annotation_file),
                annotation_hit_files_string,
                self.pathes.annotation_hit_overview_file(annotation_file)), 
             shell=True)
        call("%s %s/%s -d a %s %s > %s" % (
                self.pathes.python_bin, self.pathes.bin_folder,
                "build_annotation_table_with_read_countings.py",
                self.pathes.annotation_file(annotation_file),
                annotation_hit_files_string,
                self.pathes.annotation_hit_overview_antisense_file(annotation_file)),
             shell=True)

    def _build_annotation_hit_overview_read_normalized(self, annotation_file):
        """Create annotation hit overview table normalized by mapped
           reads.

        Arguments:
        - `annotation_file`: an (NCBI) annotation file

        """
        genome_file = self.annotation_files[annotation_file]
        mapped_reads_counting_string = ":".join(
            [str(self.grbuilder._count_mapped_reads(read_file, genome_file)) 
             for read_file in self.pathes.read_files])
        annotation_hit_files_string = " ".join(
            [self.pathes.annotation_hit_file(read_file, annotation_file) 
             for read_file in self.pathes.read_files])
        call("%s %s/%s -n %s %s %s > %s" % (
                self.pathes.python_bin, self.pathes.bin_folder,
                "build_annotation_table_with_read_countings.py",
                mapped_reads_counting_string,
                self.pathes.annotation_file(annotation_file),
                annotation_hit_files_string,
                self.pathes.annotation_hit_overview_read_normalized_file(annotation_file)), 
             shell=True)
        call("%s %s/%s -d a -n %s %s %s > %s" % (
                self.pathes.python_bin, self.pathes.bin_folder,
                "build_annotation_table_with_read_countings.py",
                mapped_reads_counting_string,
                self.pathes.annotation_file(annotation_file),
                annotation_hit_files_string,
                self.pathes._annotation_hit_overview_read_normalized_antisense_file_path(annotation_file)), 
             shell=True)

    def _build_annotation_hit_overview_nucl_normalized(self, annotation_file):
        """Create annotation hit overview table normalized by mapped
           nucleotides.

        Arguments:
        - `annotation_file`: an (NCBI) annotation file

        """
        genome_file = self.annotation_files[annotation_file]
        mapped_nucl_counting_string = ":".join(
            [str(self.grbuilder._count_mapped_nucleotides(read_file, genome_file)) 
             for read_file in self.pathes.read_files])
        annotation_hit_files_string = " ".join(
            [self.pathes.annotation_hit_file(read_file, annotation_file) 
             for read_file in self.pathes.read_files])
        call("%s %s/%s -n %s %s %s > %s" % (
                self.pathes.python_bin, self.pathes.bin_folder,
                "build_annotation_table_with_read_countings.py",
                mapped_nucl_counting_string,
                self.pathes.annotation_file(annotation_file),
                annotation_hit_files_string,
                self.pathes.annotation_hit_overview_nucl_normalized_file(annotation_file)), 
             shell=True)
        call("%s %s/%s -d a -n %s %s %s > %s" % (
                self.pathes.python_bin, self.pathes.bin_folder,
                "build_annotation_table_with_read_countings.py",
                mapped_nucl_counting_string,
                self.pathes.annotation_file(annotation_file),
                annotation_hit_files_string,
                self.pathes.annotation_hit_overview_nucl_normalized_antisense_file(annotation_file)), 
             shell=True)
