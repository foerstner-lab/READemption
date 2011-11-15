import os
import sys
sys.path.append(".")
from libs.fasta import FastaParser
from libs.grcreator import GRCreator
from libs.parameters import Parameters
from libs.paths import Paths
from libs.projectcreator import ProjectCreator
from libs.readclipper import ReadClipper
from libs.readmapper import ReadMapper
from libs.readmapperstats import ReadMapperStats
from libs.seqsizefilter import SeqSizeFilter
from libs.annotationoverlap import AnnotationOverlap

class Controller(object):

    def __init__(self, args):
        """Create an instance."""
        self.args = args
        self.paths = Paths(args.project_path)
        self.parameters = Parameters()

    def start_project(self):
        """Create a new project."""
        project_creator = ProjectCreator()
        project_creator.create_root_folder(self.args.project_path)
        project_creator.create_subfolders(self.paths.required_folders())
        project_creator.create_config_file(self.paths.config_file)
        sys.stdout.write("Created folder \"%s\" and required subfolders.\n" % (
                self.args.project_path))
        sys.stdout.write("Please copy read files into folder \"%s\" and "
                         "genome files into folder \"%s\".\n" % (
                self.paths.read_fasta_folder, self.paths.genome_folder))

    def map_reads(self):
        """Perform the mapping of the reads."""
        read_file_names = self.paths._get_read_file_names()
        genome_file_names = self.paths._get_genome_file_names()
        self.paths.set_read_files_dep_file_lists(
            read_file_names, self.parameters.min_seq_length)
        self.paths.set_genome_paths(genome_file_names)
        # TODO
        # - generate input stats before mapping
        # - read tracing after mapping
        self._prepare_reads()
        self._map_reads()
        self._generate_read_mapping_stats(read_file_names)

    def _prepare_reads(self):
        read_clipper = ReadClipper()
        read_clipper.clip(
            self.paths.read_file_paths, self.paths.clipped_read_file_paths)
        seq_size_filter = SeqSizeFilter()
        seq_size_filter.filter(
            self.paths.clipped_read_file_paths, 
            self.paths.clipped_read_file_long_enough_paths,
            self.paths.clipped_read_file_too_short_paths, 
            self.args.min_read_length)
        
    def _map_reads(self):
        read_mapper = ReadMapper(segemehl_bin=self.args.segemehl_bin)
        read_mapper.build_index(
            self.paths.genome_file_paths, self.paths.index_file_path)
        read_mapper.run_mappings(
            self.paths.clipped_read_file_long_enough_paths,
            self.paths.genome_file_paths, self.paths.index_file_path,
            self.paths.read_mapping_result_paths, 
            self.paths.unmapped_reads_paths, 
            int(self.args.threads),
            int(self.args.segemehl_accuracy),
            int(self.args.segemehl_evalue))

    def _generate_read_mapping_stats(self, read_file_names):
        ref_ids_to_file_name = self._ref_ids_to_file_name(
            self.paths.genome_file_paths)
        read_mapper_stats = ReadMapperStats()
        read_mapper_stats.count_raw_reads(
            read_file_names, self.paths.read_file_paths)
        read_mapper_stats.count_long_enough_clipped_reads(
            read_file_names, self.paths.clipped_read_file_long_enough_paths)
        read_mapper_stats.count_too_small_clipped_reads(
            read_file_names, self.paths.clipped_read_file_too_short_paths)
        read_mapper_stats.count_mappings(
            read_file_names, self.paths.read_mapping_result_paths)
        read_mapper_stats.count_unmapped_reads(
            read_file_names, self.paths.unmapped_reads_paths)
        read_mapper_stats.write_stats_to_file(
            read_file_names, ref_ids_to_file_name, 
            self.paths.read_mapping_stat_file)

    def _ref_ids_to_file_name(self, genome_file_paths):
        ref_ids_to_file_name = {}
        fasta_parser = FastaParser()
        for genome_file_path in genome_file_paths:
            genome_file = os.path.basename(genome_file_path)
            ref_seq_id = fasta_parser.header_id(
                fasta_parser.single_entry_file_header(open(genome_file_path)))
            ref_ids_to_file_name[ref_seq_id] = genome_file
        return(ref_ids_to_file_name)
    
    def create_gr_files(self):
        """Create GR files based on the combined Segemehl mappings."""
        read_file_names = self.paths._get_read_file_names()
        genome_file_names = self.paths._get_genome_file_names()
        self.paths.set_read_files_dep_file_lists(
            read_file_names, self.parameters.min_seq_length)
        self.paths.set_genome_paths(genome_file_names)
        ref_ids_to_file_name = self._ref_ids_to_file_name(
            self.paths.genome_file_paths)
        gr_creator = GRCreator()
        gr_creator.create_gr_files(
            read_file_names, self.paths.read_mapping_result_paths, 
            ref_ids_to_file_name, self.paths.gr_folder)
        gr_creator.create_read_normalized_gr_files(
            read_file_names, self.paths.read_mapping_result_paths, 
            ref_ids_to_file_name, self.paths.gr_folder_read_normalized)

    def search_annotation_overlaps(self):
        """Search for overlaps of reads and annotations."""
        read_file_names = self.paths._get_read_file_names()
        genome_file_names = self.paths._get_genome_file_names()
        self.paths.set_read_files_dep_file_lists(
            read_file_names, self.parameters.min_seq_length)
        annotation_file_names = self.paths._get_annotation_file_names()
        self.paths.set_annotation_paths(annotation_file_names)
        annotation_overlaps = AnnotationOverlap()
        annotation_overlaps.read_annotation_files(
            self.paths.annotation_file_paths)
        read_file_names = self.paths._get_read_file_names()
        annotation_overlaps.search_overlaps(
            self.paths.read_mapping_result_paths, 
            self.paths.annotation_overlap_result_paths)
        
    # def generate_report(self):
    #     """Create final report of the analysis."""
    #     self._in_project_folder()
    #     rapl_reporter = Reporter(self)
    #     report_fh = open(self.paths.report_tex_file, "w")
    #     report_fh.write(rapl_reporter.report())
    #     report_fh.close()
        
    # def _in_project_folder(self):
    #     """Check if the current directory is a RAPL project folder."""
    #     if not (os.path.exists(self.paths.config_file) and 
    #         os.path.exists(self.paths.input_folder) and 
    #         os.path.exists(self.paths.output_folder)):
    #         sys.stderr.write("Seems like the current folder is not a RAPL "
    #                          "project folder.\n")
    #         sys.stderr.write("Your are currently in \"%s\".\n" % (os.getcwd()))
    #         sys.exit(2)        

    # def _get_read_file_names(self):
    #     """Read the names of the read files."""
    #     self.read_files = sorted(os.listdir(self.paths.rna_seq_folder))
