import os
import shutil
import sys

class Rapl(object):

    def __init__(self):
        """
        
        Arguments:
        - `self`:
        """
        self._set_folder_names()

    def start_project(self, args):
        """Creates a new project
        
        Arguments:
        - `self`:
        - `args.project_name`: Name of the project root folder
        """
        self._create_root_folder(args.project_name)
        self._create_subfolders(args.project_name)
        sys.stdout.write("Created folder \"%s\" and required subfolders.\n" % (
                args.project_name)) 

    def _create_root_folder(self, project_name):
        """Create the root folder of a new project
        
        Arguments:
        - `self`:
        - `project_name`: Name of the project root folder
        """
        if not os.path.exists(project_name):
            os.mkdir(project_name)
        else:
            sys.stderr.write("Cannot create folder \"%s\"! File/folder with "
                             "the same name exists already.\n" % project_name)
            sys.exit(2)

    def _create_subfolders(self, project_name):
        """Create required subfolders
        
        Arguments:
        - `self`:
        - `project_name`: Name of the project root folder
        """
        for folder in [
            self.input_folder, self.output_folder, self.rna_seq_folder,
            self.annotation_folder, self.read_mapping_folder,
            self.read_mapping_after_clipping_folder, self.gr_folder,
            self.read_mapping_index_folder, 
            self.umapped_reads_of_first_mapping_folder,
            self.combined_mapping_folder, self.combined_mapping_split_folder,
            self.annotation_hit_folder, self.annotation_hit_overview_folder,
            self.mapping_stat_folder, self.read_tracing_folder]:
            folder_in_root_folder = "%s/%s" % (project_name, folder)
            if not os.path.exists(folder_in_root_folder):
                os.mkdir(folder_in_root_folder)

    def _set_folder_names(self):
        """Set the name of folders used in a project
        
        Arguments:
        - `self`:
        """
        self.input_folder = "input"
        self.output_folder = "output"
        self.rna_seq_folder = "%s/RNA_seqs" % self.input_folder
        self.genome_folder = "%s/genomes" % self.input_folder
        self.annotation_folder = "%s/annotation_files" % self.input_folder
        self.read_mapping_folder = "%s/read_mappings" % self.output_folder
        self.read_mapping_after_clipping_folder = (
            "%s/mappings_of_unmapped_clipped_reads" % self.output_folder)
        self.gr_folder = "%s/gr_files" % self.output_folder
        self.read_mapping_index_folder = "%s/read_mapping_index" % (
            self.output_folder)
        self.umapped_reads_of_first_mapping_folder = (
            "%s/unmapped_reads_of_first_mapping" % self.output_folder)
        self.umapped_reads_of_second_mapping_folder = (
            "%s/unmapped_reads_of_second_mapping" % self.output_folder)
        self.combined_mapping_folder = "%s/combined_read_mappings" % (
            self.output_folder)
        self.combined_mapping_split_folder = (
            "%s/combined_read_mappings_split_by_genome_files" % 
            self.output_folder)
        self.annotation_hit_folder = "%s/annotation_hits" % self.output_folder
        self.annotation_hit_overview_folder = (
            "%s/annotation_hit_overviews" % self.output_folder)
        self.mapping_stat_folder = "%s/read_mapping_stats" % (
            self.output_folder)
        self.read_tracing_folder = "%s/read_tracing" % (self.output_folder)
