import os
import shutil
import sys
import hashlib
from subprocess import call
from rapl.segemehl import SegemehlParser
from rapl.segemehl import SegemehlBuilder
from rapl.reporter import Reporter
from rapl.creator import Creator
from rapl.inputstats import InputStats
from rapl.fasta import FastaParser
from rapl.readmapper import ReadMapper
from rapl.readtracer import ReadTracer
from rapl.grbuilder import GrBuilder
from rapl.readmappingsummary import ReadMappingSummary
from rapl.annotations import Annotations
from rapl.pathes import Pathes

class Rapl(object):

    def __init__(self):
        """Create an instance."""
        self.pathes = Pathes()

    def start_project(self, args):
        """Create a new project.
        
        Arguments:
        - `args.project_name`: Name of the project root folder

        """
        rapl_creator = Creator()
        rapl_creator.create_root_folder(args.project_name)
        rapl_creator.create_subfolders(args.project_name)
        rapl_creator.create_config_file(args.project_name)
        sys.stdout.write("Created folder \"%s\" and required subfolders.\n" % (
                args.project_name))
        sys.stdout.write("Please copy read files into folder \"%s\" and "
                         "genome files into folder \"%s\".\n" % (
                self.pathes.rna_seq_folder, self.genome_folder))

    def map_reads(self):
        """Perform the mapping of the reads.

        The mapping is done using the program segemehl and takes place
        in two steps.

        """
        self._in_project_folder()
        self._get_genome_file_names()
        self._get_read_file_names()
        input_file_stats = InputStats()
        input_file_stats.create_read_file_stats()
        input_file_stats.create_genome_file_stats()
        read_mapper = ReadMapper()
        read_mapper.build_segmehl_index()
        read_mapper.run_mapping_with_raw_reads()
        read_mapper.extract_unmapped_reads_raw_read_mapping()
        read_mapper.clip_unmapped_reads()
        read_mapper.filter_clipped_reads_by_size()
        read_mapper.run_mapping_with_clipped_reads()
        read_mapper.extract_unmapped_reads_of_second_mapping()
        read_mapper.combine_mappings()
        read_mapper.filter_combined_mappings_by_a_content()
        read_mapper.split_mappings_by_genome_files()
        read_mapping_summary = ReadMappingSummary()
        read_mapping_summary.create()
        read_tracer = ReadTracer()
        read_tracer.trace_reads()
        read_tracer.create_tracing_summay()
    
    def create_gr_files(self):
        """Create GR files based on the combined Segemehl mappings. """
        self._in_project_folder()
        self._get_genome_file_names()
        self._get_read_file_names()
        gr_builder = GrBuilder()
        gr_builder.build_gr_files()
        gr_builder.build_read_normalized_gr_files()
        gr_builder.build_nucl_normalized_gr_files()

    def search_annotation_overlaps(self):
        """Search for overlaps of reads and annotations."""
        self._in_project_folder()
        self._get_genome_file_names()
        self._get_read_file_names()
        annotations = Annotations()
        annotations.find_annotation_hits()
        annotations.build_annotation_hit_overview()
        annotations.build_annotation_hit_overview_read_normalized()
        annotations.build_annotation_hit_overview_nucl_normalized()

    def generate_report(self):
        """Create final report of the analysis."""
        self._in_project_folder()
        self._get_genome_file_names()
        self._get_read_file_names()
        rapl_reporter = Reporter(self)
        report_fh = open(self.pathes.report_tex_file, "w")
        report_fh.write(rapl_reporter.report())
        report_fh.close()
        
    def _in_project_folder(self):
        """Check if the current directory is a RAPL project folder."""
        if not (os.path.exists(self.pathes.config_file) and 
            os.path.exists(self.pathes.input_folder) and 
            os.path.exists(self.pathes.output_folder)):
            sys.stderr.write("Seems like the current folder is not a RAPL "
                             "project folder.\n")
            sys.stderr.write("Your are currently in \"%s\".\n" % (os.getcwd()))
            sys.exit(2)        

    def _get_read_file_names(self):
        """Read the names of the read files."""
        self.read_files = sorted(os.listdir(self.pathes.rna_seq_folder))

    # TODO: remove
    def _set_segemehl_parameters(self):
        """Set paremeters for Segemehl."""
        self.segemehl_accuracy = 85
        self.segemehl_hit_strategy = "2"
        self.segemehl_max_e_value = 5.0
        self.segemehl_number_of_threads = 1

    # TODO: remove
    def _set_filtering_parameters(self):
        """Set parameters for sequence and Segemehl hit filtering."""
        self.min_seq_length = 12
        self.max_a_content = 70.0
        self.min_overlap = 1

    def _get_genome_file_names(self):
        """Read the names of genome files."""
        self.genome_files = sorted(os.listdir(self.pathes.genome_folder))
        
