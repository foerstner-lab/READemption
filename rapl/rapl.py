import configparser
import os
import shutil
import sys
from subprocess import call

class Rapl(object):

    def __init__(self):
        """
        
        Arguments:
        - `self`:
        """
        self._set_folder_names()
        self._set_file_names()
        self._set_bin_pathes()
        self._set_segemehl_parameters()

    def start_project(self, args):
        """Creates a new project
        
        Arguments:
        - `self`:
        - `args.project_name`: Name of the project root folder
        """
        self._create_root_folder(args.project_name)
        self._create_subfolders(args.project_name)
        self._create_config_file(args.project_name)
        sys.stdout.write("Created folder \"%s\" and required subfolders.\n" % (
                args.project_name))
        sys.stdout.write("Please copy read files into folder \"%s\" and "
                         "genome files into folder \"%s\".\n" % (
                self.rna_seq_folder, self.genome_folder))

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
            self.read_mapping_index_folder, self.genome_folder,
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

    def _set_file_names(self):
        """Set name of common files"""
        self.config_file = "rapl.config"

    def _set_bin_pathes(self):
        self.segemehl_bin = "segemehl"

    def _set_segemehl_parameters(self):
        self.segemehl_accuracy = 85
        self.segemehl_hit_strategy = "1"
        self.segemehl_max_e_value = 10
        self.segemehl_number_of_threads = 1
        
    def _create_config_file(self, project_name):
        """Creates a config file
        
        Arguments:
        - `self`:
        - `project_name`: Name of the project root folder
        """
        config = configparser.RawConfigParser()
        config.add_section('RAPL')
        with open("%s/%s" % (project_name, self.config_file), 
                  "w") as configfile:
            config.write(configfile)

    def map_reads(self, args):
        """Perform the mapping of the reads

        The mapping is done using the program segemehl and takes place
        in two steps.

        Arguments:
        - `self`:
        - `args`: 
        """
        self._in_project_folder()
        self._get_genome_file_names()
        self._get_read_file_names()
        self._build_segmehl_index()
        self._run_mapping_with_raw_reads()
        sefl._extract_unmapped_reads_of_second_mapping()

    def _in_project_folder(self):
        """Check if the current directory is a RAPL project folder"""
        if not (os.path.exists(self.config_file) and 
            os.path.exists(self.input_folder) and 
            os.path.exists(self.output_folder)):
            sys.stderr.write("Seems like the current folder is not a RAPL "
                             "project folder.\n")
            sys.stderr.write("Your are currently in \"%s\".\n" % (os.getcwd()))
            sys.exit(2)        

    def _get_read_file_names(self):
        """Read the name of the read files"""
        self.read_files = os.listdir(self.rna_seq_folder)

    def _get_genome_file_names(self):
        """Read the names of genome files"""
        self.genome_files = os.listdir(self.genome_folder)
        
    def _build_segmehl_index(self):
        """Create the segemehl index based on the genome files"""
        call("echo %s -x %s -d %s" % (
                self.segemehl_bin, self._segemehl_index_path(),
                " ".join(self._genome_file_pathes())), 
             shell=True)

    def _segemehl_index_path(self):
        """Full path the the in segemehl index file"""
        return("%s/%s"  % (
                self.read_mapping_index_folder, self._segemehl_index_name()))

    def _genome_file_path(self, genome_file):
        """Full path of a given genome file

        Arguments:
        - `self`:
        - `genome_file`: genome file name
        """
        return("%s/%s" % (self.genome_folder, genome_file))

    def _segemehl_index_name(self):
        """Name of the segemehl index file"""
        index_file_name = "_".join(self.genome_files) + ".idx"
        index_file_name.replace(".fa", "")
        return(index_file_name)

    def _genome_file_pathes(self):
        """Full pathes of all genome files"""
        return([self._genome_file_path(genome_file) 
                for genome_file in self.genome_files])
    
    def _run_mapping_with_raw_reads(self):
        """Run the mapping of the raw reads."""
        for read_file in self.read_files:
            self._run_segemehl_search(
                self._read_file_path(read_file),
                self._raw_read_mapping_output_path(read_file),
                self._unmapped_raw_read_file_path(read_file))

    def _raw_read_mapping_output_path(self, read_file):
        """Full path of the output file of a segemehl run

        Arguments:
        - `self`:
        - `read_file`: read file name that is mapped
        """
        return("%s/%s_mapped_to_%s" % (
                self.read_mapping_folder, read_file, self._segemehl_index_name()))

    def _run_segemehl_search(self, read_file_path, output_file_path, 
                             unmapped_read_file_path):
        """Call segemehl to do a mapping

        Arguments:
        - `self`:
        - `read_file_path`:
        - `output_file_path`:
        - `unmapped_read_file_path`:
        """
        call("%s -E %s -H %s -A %s -t %s -i %s -d %s -q %s -o %s" % (
                self.segemehl_bin,
                self.segemehl_max_e_value,
                self.segemehl_hit_strategy,
                self.segemehl_accuracy,
                self.segemehl_number_of_threads,
                self._segemehl_index_path(),
                " ".join(self._genome_file_pathes()),
                read_file_path,
                output_file_path),
             shell=True)
        # FOR THE UPCOMING VERSION
        # call("%s -E %s -H %s -A %s -t %s -i %s -d %s -q %s -o %s -u %s" % (
        #         self.segemehl_bin,
        #         self.segemehl_max_e_value,
        #         self.segemehl_hit_strategy,
        #         self.segemehl_accuracy,
        #         self.segemehl_number_of_threads,
        #         self._segemehl_index_path(),
        #         self._genome_file_pathes(),
        #         read_file_path,
        #         output_file_path,
        #         unmapped_read_file_path),
        #      shell=True)

    def _read_file_path(self, read_file):
        """The full path of a given read file

        Arguments:
        - `self`:
        - `read_file`:
        """
        return("%s/%s" % (self.rna_seq_folder, read_file))

    def _unmapped_raw_read_file_path(self, read_file):
        return("%s/%s.unmapped.fa" % (
                self.umapped_reads_of_first_mapping_folder, read_file))
