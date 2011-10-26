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

class Controller(object):

    def __init__(self):
        """Create an instance."""
        self.paths = Paths()
        self.parameters = Parameters()

    def start_project(self, args):
        """Create a new project.
        
        Arguments:
        - `args.project_name`: Name of the project root folder

        """
        project_creator = ProjectCreator()
        project_creator.create_root_folder(args.project_name)
        project_creator.create_subfolders(
            args.project_name, self.paths.required_folders())
        project_creator.create_config_file(
            args.project_name, self.paths.config_file)
        sys.stdout.write("Created folder \"%s\" and required subfolders.\n" % (
                args.project_name))
        sys.stdout.write("Please copy read files into folder \"%s\" and "
                         "genome files into folder \"%s\".\n" % (
                self.paths.read_fasta_folder, self.paths.genome_folder))

    def map_reads(self, args=None):
        """Perform the mapping of the reads."""
        read_file_names = self.paths._get_read_file_names()
        genome_file_names = self.paths._get_genome_file_names()
        self.paths.set_read_files_dep_file_lists(
            read_file_names, self.parameters.min_seq_length)
        self.paths.set_genome_paths(genome_file_names)
        ref_ids_to_file_name = self._ref_ids_to_file_name(
            self.paths.genome_file_paths)
        # self._in_project_folder()
        # input_file_stats = InputStats()
        # input_file_stats.create_read_file_stats()
        # input_file_stats.create_genome_file_stats()
        read_clipper = ReadClipper()
        read_clipper.clip(
            self.paths.read_file_paths, self.paths.clipped_read_file_paths)
        seq_size_filter = SeqSizeFilter()
        seq_size_filter.filter(
            self.paths.clipped_read_file_paths, 
            self.paths.clipped_read_file_long_enough_paths,
            self.paths.clipped_read_file_too_short_paths, 
            args.min_read_length)
        read_mapper = ReadMapper(segemehl_bin=args.segemehl_bin)
        read_mapper.build_index(
            self.paths.genome_file_paths, self.paths.index_file_path)
        read_mapper.run_mappings(
            self.paths.clipped_read_file_long_enough_paths,
            self.paths.genome_file_paths, self.paths.index_file_path,
            self.paths.read_mapping_result_paths, 
            self.paths.unmapped_reads_paths, 
            int(args.threads),
            int(args.segemehl_accuracy),
            int(args.segemehl_evalue))
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

        # read_mapper.select_uniquely_mapped_reads()
        # read_tracer = ReadTracer()
        # read_tracer.trace_reads()
        # read_tracer.create_tracing_summay()
        # read_tracer_viz = ReadTracerViz()
        # read_tracer_viz.create_mapping_length_histograms()
    
    def create_gr_files(self):
        """Create GR files based on the combined Segemehl mappings."""
        #self._in_project_folder()
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

    # def search_annotation_overlaps(self):
    #     """Search for overlaps of reads and annotations."""
    #     self._in_project_folder()
    #     annotations = Annotations()
    #     annotations.find_annotation_hits()
    #     annotations.quantify_mapping_redundancy()
    #     annotations.build_annotation_hit_overview()
    #     annotations.build_annotation_hit_overview_read_normalized()
    #     annotations.build_annotation_hit_overview_nucleotide_normalized()
    #     annotations.build_annotation_hit_overview_rpkm_normalized()
    #     annotations.count_reads_in_intergenic_regions()
        
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
