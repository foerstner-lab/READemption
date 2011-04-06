import os
import shutil
import sys
import hashlib
from subprocess import call
from rapl.segemehl import SegemehlParser
from rapl.segemehl import SegemehlBuilder
from rapl.raplreporter import RaplReporter
from rapl.raplcreator import RaplCreator
from rapl.raplinputstats import RaplInputStats
from rapl.fasta import FastaParser
from rapl.raplreadmapper import RaplReadMapper
from rapl.raplreadtracer import RaplReadTracer
from rapl.raplgrbuilder import RaplGrBuilder
from rapl.raplreadmappingsummary import ReadMappingSummary

class Rapl(object):

    def __init__(self):
        """Create an instance."""
        self._set_folder_names()
        self._set_file_names()
        self._set_segemehl_parameters()
        self._set_filtering_parameters()

    def start_project(self, args):
        """Create a new project.
        
        Arguments:
        - `args.project_name`: Name of the project root folder

        """
        rapl_creator = RaplCreator()
        rapl_creator.create_root_folder(args.project_name)
        rapl_creator.create_subfolders(args.project_name)
        rapl_creator.create_config_file(args.project_name)
        sys.stdout.write("Created folder \"%s\" and required subfolders.\n" % (
                args.project_name))
        sys.stdout.write("Please copy read files into folder \"%s\" and "
                         "genome files into folder \"%s\".\n" % (
                self.rna_seq_folder, self.genome_folder))

    def map_reads(self, args):
        """Perform the mapping of the reads.

        The mapping is done using the program segemehl and takes place
        in two steps.

        Arguments:
        - `args`: command line arguments
        
        """
        self._in_project_folder()
        self._get_genome_file_names()
        self._get_read_file_names()
        input_file_stats = RaplInputStats()
        input_file_stats.create_read_file_stats()
        input_file_stats.create_genome_file_stats()
        read_mapper = RaplReadMapper()
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
    
    def create_gr_files(self, args):
        """Create GR files based on the combined Segemehl mappings.

        Arguments:
        - `args`: command line arguments

        """
        self._in_project_folder()
        self._get_genome_file_names()
        self._get_read_file_names()
        gr_builder = RaplGrBuilder()
        gr_builder.build_gr_files()
        gr_builder.build_read_normalized_gr_files()
        gr_builder.build_nucl_normalized_gr_files()

    def search_annotation_overlaps(self, args):
        """Search for overlaps of reads and annotations.

        Arguments:
        - `args`: command line arguments

        """
        self._read_config_file()
        self._in_project_folder()
        self._get_genome_file_names()
        self._get_read_file_names()
        self._get_annotation_files_from_config()
        self.find_annotation_hits()
        self.build_annotation_hit_overview()
        self.build_annotation_hit_overview_read_normalized()
        self.build_annotation_hit_overview_nucl_normalized()

    def generate_report(self, args):
        """Create final report of the analysis.

        Arguments:
        - `args`: command line arguments

        """
        self._read_config_file()
        self._in_project_folder()
        self._get_genome_file_names()
        self._get_read_file_names()
        self._get_annotation_files_from_config()
        rapl_reporter = RaplReporter(self)
        report_fh = open(self.report_tex_file, "w")
        report_fh.write(rapl_reporter.report())
        report_fh.close()
        
    def _in_project_folder(self):
        """Check if the current directory is a RAPL project folder."""
        if not (os.path.exists(self.config_file) and 
            os.path.exists(self.input_folder) and 
            os.path.exists(self.output_folder)):
            sys.stderr.write("Seems like the current folder is not a RAPL "
                             "project folder.\n")
            sys.stderr.write("Your are currently in \"%s\".\n" % (os.getcwd()))
            sys.exit(2)        

    def _get_read_file_names(self):
        """Read the names of the read files."""
        self.read_files = sorted(os.listdir(self.rna_seq_folder))

    def _find_annotation_hits(self, read_file, annotation_file):
        """Search for overlaps of reads and annotations.

        Arguments:
        - `read_file`: orignal read file used to generate the mappings.
        - `annotation_file`: an (NCBI) annotation file

        """
        genome_file = self.annotation_files[annotation_file]
        call("%s %s/segemehl_hit_annotation_mapping.py -m %s -o %s %s %s" % (
                self.python_bin, self.bin_folder,
                self.min_overlap,
                self._annotation_hit_file_path(read_file, annotation_file),
                self._combined_mapping_file_a_filtered_split_path(
                    read_file, genome_file),
                self._annotation_file_path(annotation_file)), shell=True)

    def _get_annotation_files_from_config(self):
        """Get the annations files from the config files.

        It extracts a dictionary that contains the names of the
        annotation files as keys and the names of the corresponding
        genome files as values.

        """
        self.annotation_files = self.config["annotation_and_genomes_files"]

    def _read_config_file(self):
        """Read the config file."""
        self.config = json.loads(open(self.config_file).read())

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
            [self._annotation_hit_file_path(read_file, annotation_file) 
             for read_file in self.read_files])
        call("%s %s/%s %s %s > %s" % (
                self.python_bin, self.bin_folder,
                "build_annotation_table_with_read_countings.py",
                self._annotation_file_path(annotation_file),
                annotation_hit_files_string,
                self._annotation_hit_overview_file_path(annotation_file)), 
             shell=True)
        call("%s %s/%s -d a %s %s > %s" % (
                self.python_bin, self.bin_folder,
                "build_annotation_table_with_read_countings.py",
                self._annotation_file_path(annotation_file),
                annotation_hit_files_string,
                self._annotation_hit_overview_antisense_file_path(annotation_file)),
             shell=True)

    def _build_annotation_hit_overview_read_normalized(self, annotation_file):
        """Create annotation hit overview table normalized by mapped
           reads.

        Arguments:
        - `annotation_file`: an (NCBI) annotation file

        """
        genome_file = self.annotation_files[annotation_file]
        mapped_reads_counting_string = ":".join(
            [str(self._count_mapped_reads(read_file, genome_file)) 
             for read_file in self.read_files])
        annotation_hit_files_string = " ".join(
            [self._annotation_hit_file_path(read_file, annotation_file) 
             for read_file in self.read_files])
        call("%s %s/%s -n %s %s %s > %s" % (
                self.python_bin, self.bin_folder,
                "build_annotation_table_with_read_countings.py",
                mapped_reads_counting_string,
                self._annotation_file_path(annotation_file),
                annotation_hit_files_string,
                self._annotation_hit_overview_read_normalized_file_path(annotation_file)), 
             shell=True)
        call("%s %s/%s -d a -n %s %s %s > %s" % (
                self.python_bin, self.bin_folder,
                "build_annotation_table_with_read_countings.py",
                mapped_reads_counting_string,
                self._annotation_file_path(annotation_file),
                annotation_hit_files_string,
                self._annotation_hit_overview_read_normalized_file_path(annotation_file)), 
             shell=True)

    def _build_annotation_hit_overview_nucl_normalized(self, annotation_file):
        """Create annotation hit overview table normalized by mapped
           nucleotides.

        Arguments:
        - `annotation_file`: an (NCBI) annotation file

        """
        genome_file = self.annotation_files[annotation_file]
        mapped_nucl_counting_string = ":".join(
            [str(self._count_mapped_nucleotides(read_file, genome_file)) 
             for read_file in self.read_files])
        annotation_hit_files_string = " ".join(
            [self._annotation_hit_file_path(read_file, annotation_file) 
             for read_file in self.read_files])
        call("%s %s/%s -n %s %s %s > %s" % (
                self.python_bin, self.bin_folder,
                "build_annotation_table_with_read_countings.py",
                mapped_nucl_counting_string,
                self._annotation_file_path(annotation_file),
                annotation_hit_files_string,
                self._annotation_hit_overview_nucl_normalized_file_path(annotation_file)), 
             shell=True)
        call("%s %s/%s -d a -n %s %s %s > %s" % (
                self.python_bin, self.bin_folder,
                "build_annotation_table_with_read_countings.py",
                mapped_nucl_counting_string,
                self._annotation_file_path(annotation_file),
                annotation_hit_files_string,
                self._annotation_hit_overview_nucl_normalized_file_path(annotation_file)), 
             shell=True)

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
    

    ####################        
    # Paths
    ####################
        
    def _read_file_path(self, read_file):
        """Return the full path of a given read file.

        Arguments:
        - `read_file`: name of the read file
        """
        return("%s/%s" % (self.rna_seq_folder, read_file))

    def _unmapped_raw_read_file_path(self, read_file):
        """Return the full path of unmapped reads of the first run.

        Arguments:
        - `read_file`: name of the read file
        """
        return("%s/%s.unmapped.fa" % (
                self.umapped_reads_of_first_mapping_folder, read_file))

    def _segemehl_index_path(self):
        """Return the full path the the in segemehl index file."""
        return("%s/%s"  % (
                self.read_mapping_index_folder, self._segemehl_index_name()))

    def _genome_file_path(self, genome_file):
        """Return the ull path of a given genome file

        Arguments:
        - `genome_file`: genome file name
        """
        return("%s/%s" % (self.genome_folder, genome_file))

    def _genome_file_paths(self):
        """Return the full paths of all genome files"""
        return([self._genome_file_path(genome_file) 
                for genome_file in self.genome_files])

    def _raw_read_mapping_output_path(self, read_file):
        """Return the full path of the output file of a segemehl run

        Arguments:
        - `read_file`: read file name that is mapped
        """
        return("%s/%s_mapped_to_%s" % (
                self.read_mappings_first_run_folder, read_file, self._segemehl_index_name()))

    def _unmapped_read_clipped_path(self, read_file):
        """Return the full path of a file with clipped reads

        Arguments:
        - `read_file`: name of the read file
        """
        return("%s/%s.unmapped.fa.clipped.fa" % (
                self.umapped_reads_of_first_mapping_folder, read_file))


    def _unmapped_clipped_size_filtered_read_path(self, read_file):
        """Return the full path of clipped and size filtered reads.

        Arguments:
        - `read_file`: name of the read file

        """
        return("%s/%s.unmapped.fa.clipped.fa.size_filtered_gtoe_%sbp.fa" % (
                self.umapped_reads_of_first_mapping_folder,
                read_file, self.min_seq_length))

    def _clipped_reads_mapping_output_path(self, read_file):
        """Return the path of the mapping file of the second run.

        Arguments:
        - `read_file`: name of the read file
        """
        return("%s/%s.clipped_mapped_to_%s" % (
                self.read_mappings_second_run_folder,
                read_file, self._segemehl_index_name()))

    def _unmapped_reads_of_clipped_reads_file_path(self, read_file):
        """Return the full path of the unmapped reads of the 2nd run.

        Arguments:
        - `read_file`: name of the read file
        """
        return("%s/%s.unmapped.fa"  % (self.umapped_reads_of_first_mapping_folder, 
                           read_file))

    def _unmapped_reads_second_mapping_path(self, read_file):
        """Return the path of unmapped reads of the second run.

        Arguments:
        - `read_file`: name of the read file
        """
        return("%s/%s.unmapped.fa" % (
                self.umapped_reads_of_second_mapping_folder, read_file))

    def _combined_mapping_file_path(self, read_file):
        """Return the path of the combined mappings.

        Arguments:
        - `read_file`: name of the read file
        """
        return("%s/%s_mapped_to_%s.combined" % (
                self.combined_mappings_folder, read_file,
                self._segemehl_index_name()))

    def _combined_mapping_file_a_filtered_split_path(self, read_file, genome_file):
        """Return the path of the split filtered combined mappings.
        
        Arguments:
        - `read_file`: name of the read file
        """
        return("%s/%s_mapped_to_%s.combined.filtered_ltoe_%s%%_A.txt.from_%s_only" % (
                self.combined_mapping_split_folder, read_file, 
                self._segemehl_index_name(), self.max_a_content, genome_file))

    def _combined_mapping_file_a_filtered_path(self, read_file):
        """Return the path of the filtered combined mappings.

        Arguments:
        - `read_file`: name of the read file
        """
        return("%s.filtered_ltoe_%s%%_A.txt" % (
                self._combined_mapping_file_path(read_file),
                self.max_a_content))

    def _unmapped_clipped_size_failed_read_path(self, read_file):
        """Return the path of size filter failed clipped reads.

        Arguments:
        - `read_file`: name of the read file
        """
        return("%s/%s.unmapped.fa.clipped.fa.size_filtered_lt_%sbp.fa" % (
                self.umapped_reads_of_first_mapping_folder,
                read_file, self.min_seq_length))

    def _combined_mapping_file_a_filter_failed_path(self, read_file):
        """Return the path of the A-content filter failed reads.

        Arguments:
        - `read_file`: name of the read file
        """
        return("%s.filtered_gt_%s%%_A.txt" % (
                self._combined_mapping_file_path(read_file),
                self.max_a_content))

    def _trace_file_path(self, read_file):
        """Return the path of the trace file of a read file.

        Arguments:
        - `read_file`: name of the read file
        """
        return("%s/%s.mapping_tracing.csv" % (
                self.read_tracing_folder, read_file))

    def _gr_file_path(self, read_file, genome_file):
        """Return the GR file path of a given read and genome file.

        Arguments:
        - `read_file,`: name of the read file
        - `genome_file`: name of the genome file
        """
        return("%s/%s_in_%s.gr" % (self.gr_folder, read_file, genome_file))

    def _gr_read_normalized_file_path(self, read_file, genome_file):
        """Return the GR read normalized file path of a given read and
        genome file.

        Arguments:
        - `read_file,`: name of the read file
        - `genome_file`: name of the genome file
        """
        return("%s/%s_in_%s.gr" % (
                self.gr_folder_read_normalized, read_file, genome_file))

    def _gr_nucl_normalized_file_path(self, read_file, genome_file):
        """Return the GR nucleotide normalized file path of a given read and
        genome file.

        Arguments:
        - `read_file,`: name of the read file
        - `genome_file`: name of the genome file
        """
        return("%s/%s_in_%s.gr" % (
                self.gr_folder_nucl_normalized, read_file, genome_file))

    def _annotation_hit_file_path(self, read_file, annotation_file):
        """Return the path of the annoation hit file.

        Arguments:
        - `read_file,`: name of the read file
        - `annotation_file`: name of the (NCBI) annotation file
        """
        return("%s/%s_in_%s_annotation_hits" % (
                self.annotation_hit_folder, read_file, annotation_file))

    def _annotation_file_path(self, annotation_file):
        """Return the path of a given annotation files.

        Arguments:
        - `annotation_file`: name of the (NCBI) annotation file
        """
        return("%s/%s" % (self.annotation_folder , annotation_file))

    def _annotation_hit_overview_file_path(self, annotation_file):
        """Return the path of the annotation overview file.

        Arguments:
        - `annotation_file`: name of the (NCBI) annotation file
        """
        return("%s/%s_all_annotation_hits_sense.csv" % (
                self.annotation_hit_overview_folder, annotation_file))

    def _annotation_hit_overview_antisense_file_path(self, annotation_file):
        """Return the path of the annotation overview file for antisense hits.

        Arguments:
        - `annotation_file`: name of the (NCBI) annotation file
        """
        return("%s/%s_all_annotation_hits_antisense.csv" % (
                self.annotation_hit_overview_folder, annotation_file))

    def _annotation_hit_overview_read_normalized_file_path(self, annotation_file):
        """Return the path of the annotation overview normalized by
           mapped reads file.
        """
        return("%s/%s_all_annotation_hits_normalized_by_reads_sense.csv" % (
                self.annotation_hit_overview_read_normalized_folder, 
                annotation_file))

    def _annotation_hit_overview_read_normalized_antisense_file_path(self, annotation_file):
        """Return the path of the annotation overview normalized by
           mapped reads file.
        """
        return("%s/%s_all_annotation_hits_normalized_by_reads_antisense.csv" % (
                self.annotation_hit_overview_read_normalized_folder, 
                annotation_file))

    def _annotation_hit_overview_nucl_normalized_file_path(self, annotation_file):
        """Return the path of the annotation overview normalized by
           mapped nucleotides file.
        """
        return("%s/%s_all_annotation_hits_normalized_by_nucleotides_sense.csv" % (
                self.annotation_hit_overview_nucl_normalized_folder, 
                annotation_file))

    def _annotation_hit_overview_nucl_normalized_antisense_file_path(self, annotation_file):
        """Return the path of the annotation overview normalized by
           mapped nucleotides file.
        """
        return("%s/%s_all_annotation_hits_normalized_by_nucleotides_antisense.csv" % (
                self.annotation_hit_overview_nucl_normalized_folder, 
                annotation_file))



    def _set_folder_names(self):
        """Set the name of folders used in a project."""
        self.input_folder = "input"
        self.output_folder = "output"
        self.rna_seq_folder = "%s/RNA_seqs" % self.input_folder
        self.genome_folder = "%s/genomes" % self.input_folder
        self.input_file_stats_folder = "%s/input_file_stats" % self.output_folder
        self.annotation_folder = "%s/annotation_files" % self.input_folder
        self.read_mappings_first_run_folder = "%s/read_mappings_first_run" % (
            self.output_folder)
        self.read_mappings_second_run_folder = (
            "%s/read_mappings_second_run" % self.output_folder)
        self.gr_folder = "%s/gr_files" % self.output_folder
        self.gr_folder_read_normalized = "%s/gr_read_normalized_files" % (
            self.output_folder)
        self.gr_folder_nucl_normalized = "%s/gr_nucleotide_normalized_files" % (
            self.output_folder)
        self.read_mapping_index_folder = "%s/read_mapping_index" % (
            self.output_folder)
        self.umapped_reads_of_first_mapping_folder = (
            "%s/unmapped_reads_of_first_mapping" % self.output_folder)
        self.umapped_reads_of_second_mapping_folder = (
            "%s/unmapped_reads_of_second_mapping" % self.output_folder)
        self.combined_mappings_folder = "%s/read_mappings_combined" % (
            self.output_folder)
        self.combined_mapping_split_folder = (
            "%s/read_mappings_combined_split_by_genome_files" % 
            self.output_folder)
        self.annotation_hit_folder = "%s/annotation_hits" % self.output_folder
        self.annotation_hit_overview_folder = (
            "%s/annotation_hit_overviews" % self.output_folder)
        self.annotation_hit_overview_read_normalized_folder = (
            "%s/annotation_hit_overviews_read_normalized" % self.output_folder)
        self.annotation_hit_overview_nucl_normalized_folder = (
            "%s/annotation_hit_overviews_nucleotide_normalized" % 
            self.output_folder)
        self.read_tracing_folder = "%s/read_tracing" % (self.output_folder)
        self.report_folder = "%s/reports_and_stats" % (self.output_folder)
        # Currently not needed
        #self.mapping_stat_folder = "%s/read_mapping_stats" % (
        #    self.output_folder)

    def _set_file_names(self):
        """Set name of common files."""
        self.config_file = "rapl.config"
        self.read_file_stats = "%s/read_file_stats.txt" % (
            self.input_file_stats_folder)
        self.genome_file_stats = "%s/genome_file_stats.txt" % (
            self.input_file_stats_folder)
        self.annotation_file_stats = "%s/annotation_file_stats.txt" % (
            self.input_file_stats_folder)
        self.tracing_summary_file = "%s/read_tracing_summary.csv" % (
            self.report_folder)
        self.report_tex_file = "%s/report.tex" % (
            self.report_folder)
        self.lib_genome_read_mapping_summary = (
            "%s/read_coutings_per_genome_file.cvs" % (self.report_folder))


    def _get_genome_file_names(self):
        """Read the names of genome files."""
        self.genome_files = sorted(os.listdir(self.genome_folder))
        
