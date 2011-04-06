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

class Rapl(object):

    def __init__(self):
        """Create an instance."""
        self._set_folder_names()
        self._set_file_names()
        self._set_bin_paths()
        self._set_segemehl_parameters()
        self._set_filtering_parameters()

    def start_project(self, args):
        """Create a new project.
        
        Arguments:
        - `args.project_name`: Name of the project root folder

        """
        # OBSOLETE - CLEAN
        #self._create_root_folder(args.project_name)
        #self._create_subfolders(args.project_name)
        #self._create_config_file(args.project_name)
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
        # self._create_read_file_stats()
        # self._create_genome_file_stats()
        
        self.build_segmehl_index()
        self.run_mapping_with_raw_reads()
        self.extract_unmapped_reads_raw_read_mapping()
        self.clip_unmapped_reads()
        self.filter_clipped_reads_by_size()
        self.run_mapping_with_clipped_reads()
        self.extract_unmapped_reads_of_second_mapping()
        self.combine_mappings()
        self.filter_combined_mappings_by_a_content()
        self.split_mappings_by_genome_files()
        self.trace_reads_after_mapping()
        self.create_tracing_summay()
        self.create_lib_genome_summary()


    def _set_bin_paths(self):
        """Set file/folder paths for some needed binaries."""
        self.segemehl_bin = "segemehl"
        # DEV
        self.python_bin = "/opt/Python-3.2/python"
        self.bin_folder = "~/rapl_tools"
        
    def _set_segemehl_parameters(self):
        """Set paremeters for Segemehl."""
        self.segemehl_accuracy = 85
        self.segemehl_hit_strategy = "2"
        self.segemehl_max_e_value = 5.0
        self.segemehl_number_of_threads = 1

    def _set_filtering_parameters(self):
        """Set parameters for sequence and Segemehl hit filtering."""
        self.min_seq_length = 12
        self.max_a_content = 70.0
        self.min_overlap = 1
        
    
    def create_gr_files(self, args):
        """Create GR files based on the combined Segemehl mappings.

        Arguments:
        - `args`: command line arguments

        """
        self._in_project_folder()
        self._get_genome_file_names()
        self._get_read_file_names()
        self.build_gr_files()
        self.build_read_normalized_gr_files()
        self.build_nucl_normalized_gr_files()

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

    def _get_genome_file_names(self):
        """Read the names of genome files."""
        self.genome_files = sorted(os.listdir(self.genome_folder))
        
    def build_segmehl_index(self):
        """Create the segemehl index based on the genome files."""
        call("%s -x %s -d %s" % (
                self.segemehl_bin, self._segemehl_index_path(),
                " ".join(self._genome_file_paths())), 
             shell=True)

    def _segemehl_index_name(self):
        """Return the name of the segemehl index file."""
        # TODO Avoid too long file name later.
        #index_file_name = "_".join(self.genome_files) + ".idx"
        #index_file_name.replace(".fa", "")
        index_file_name = "genome.idx"
        return(index_file_name)

    def run_mapping_with_raw_reads(self):
        """Run the mapping of the raw reads using segemehl"""
        for read_file in self.read_files:
            self._run_segemehl_search(
                self._read_file_path(read_file),
                self._raw_read_mapping_output_path(read_file),
                self._unmapped_raw_read_file_path(read_file))

    def _run_segemehl_search(self, read_file_path, output_file_path, 
                             unmapped_read_file_path):
        """Call segemehl to do a mapping.

        Arguments:
        - `read_file_path`: the file path of the read fasta file
        - `output_file_path`: the path of the Segemehl output
        - `unmapped_read_file_path`: the path of the fasta of 
                                     unmapped reads

        """
        call("%s -E %s -H %s -A %s -t %s -i %s -d %s -q %s -o %s" % (
                self.segemehl_bin,
                self.segemehl_max_e_value,
                self.segemehl_hit_strategy,
                self.segemehl_accuracy,
                self.segemehl_number_of_threads,
                self._segemehl_index_path(),
                " ".join(self._genome_file_paths()),
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
        #         self._genome_file_paths(),
        #         read_file_path,
        #         output_file_path,
        #         unmapped_read_file_path),
        #      shell=True)

    def extract_unmapped_reads_raw_read_mapping(self):
        """Extract unmapped reads of first mapping round."""
        for read_file in self.read_files:
            self._extract_unmapped_reads(
                self._read_file_path(read_file),
                self._raw_read_mapping_output_path(read_file),
                self._unmapped_raw_read_file_path(read_file))

    def _extract_unmapped_reads(self, read_file, mapping_file, 
                                output_read_file):
        """Extract unmapped reads of a Segemehl mapping run.

        Arguments:
        - `read_file`: name of the read fasta file
        - `mapping_file,`: name of the Segemehl mapping file
        - `output_read_file`: name of the fasta file with the
                              unmapped reads

        """
        call("%s %s/%s -o %s %s %s" % (
                self.python_bin, self.bin_folder, "extract_unmapped_fastas.py", 
                output_read_file, read_file, mapping_file), shell=True)

    def clip_unmapped_reads(self):
        """Clip reads unmapped in the first segemehl run."""
        for read_file in self.read_files:
            self._clip_reads(self._unmapped_raw_read_file_path(read_file))

    def _clip_reads(self, unmapped_raw_read_file_path):
        """Remove the poly-A tail of reads in a file.

        Arguments:
        - `unmapped_raw_read_file_path`: path of the fasta file that
                                         contains unmapped reads

        """
        call("%s %s/poly_a_clipper.py %s" % (self.python_bin,
                self.bin_folder, unmapped_raw_read_file_path), shell=True)

    def filter_clipped_reads_by_size(self):
        """Filter clipped reads sequence length.

        For each read file two output files are generated. One
        contains reads with a size equal or higher than the given
        cut-off. One for the smaller ones.

        """
        for read_file in self.read_files:
            self._filter_reads_by_size(self._unmapped_read_clipped_path(read_file))

    def _filter_reads_by_size(self, read_file_path):
        """Filter reads by sequence length.

        Arguments:
        - `read_file_path`: path of the fasta file that will be split.

        """
        call("%s %s/filter_fasta_entries_by_size.py %s %s" % (
                self.python_bin, self.bin_folder, read_file_path, 
                self.min_seq_length), shell=True)

    def run_mapping_with_clipped_reads(self):
        """Run the mapping with clipped and size filtered reads."""
        for read_file in self.read_files:
            self._run_segemehl_search(
                self._unmapped_clipped_size_filtered_read_path(read_file), 
                self._clipped_reads_mapping_output_path(read_file),
                self._unmapped_reads_of_clipped_reads_file_path(read_file))

    def extract_unmapped_reads_of_second_mapping(self):
        """Extract reads that are not mapped in the second run."""
        for read_file in self.read_files:
            self._extract_unmapped_reads(
                self._unmapped_clipped_size_filtered_read_path(read_file),
                self._clipped_reads_mapping_output_path(read_file),
                self._unmapped_reads_second_mapping_path(read_file))

    def combine_mappings(self):
        """Combine the results of both segemehl mappings for all libraries."""
        for read_file in self.read_files:
            self._combine_mappings(read_file)

    def _combine_mappings(self, read_file):
        """Combine the results of both segemehl mappings.

        Arguments:
        - `read_file`: the name of the read file that was used to generate
                       the Segemehl mappings.

        """
        comined_mappings_fh = open(self._combined_mapping_file_path(read_file), "w")
        comined_mappings_fh.write(open(self._raw_read_mapping_output_path(read_file)).read())
        comined_mappings_fh.write(open(self._clipped_reads_mapping_output_path(read_file)).read())
        comined_mappings_fh.close()

    def filter_combined_mappings_by_a_content(self):
        """Filter Segemehl mapping file entries by amount of A content.

        This removes sequences that exceed a certain amount of A that
        might be introduced during the sample preparion process.

        """
        for read_file in self.read_files:
            self._filter_combined_mappings_by_a_content(read_file)
    
    def _filter_combined_mappings_by_a_content(self, mapping_file):
        """Filter Segemehl mapping file entries by A-content.

        Two files are produced. One that contains reads that have an
        A-content higher than the cut-off value, one that contains the
        reads that have A-content equal or lower than the cut-off
        value.

        Arguments:
        - `mapping_file`: the input mapping file 

        """
        call("%s %s/filter_segemehl_by_nucleotide_percentage.py %s A %s " % (
            self.python_bin, self.bin_folder, 
            self._combined_mapping_file_path(mapping_file), self.max_a_content), 
             shell=True)

    def split_mappings_by_genome_files(self):
        """Split the Segemehl result entries by genome file."""
        headers_of_genome_files = self._headers_of_genome_files()
        for read_file in self.read_files:
            self._split_mapping_by_genome_files(
                read_file, headers_of_genome_files)

    def _split_mapping_by_genome_files(self, read_file, headers_of_genome_files):
        """Split the Segemehl results by the target genome files.

        Arguments:
        - `read_file,`: the read file that was used to generate the combined
                        Segemehl mapping file
        - `headers_of_genome_files`: A dictionary that contains the headers
                                     of the genome files as keys and the
                                     name of their files as values.

        """
        segemehl_parser = SegemehlParser()
        segemehl_builder = SegemehlBuilder()
        file_handles = {}
        # Open an output file for each genome file. Needed as some
        # genome files don't have any mapping and so their mapping
        # file would not be created otherwise and be missing later.
        for genome_file in self.genome_files:
            output_file = self._combined_mapping_file_a_filtered_split_path(
                read_file, genome_file)
            file_handles["%s-%s" % (read_file, genome_file)] = open(
                output_file, "w")
        for entry in segemehl_parser.entries(
            self._combined_mapping_file_a_filtered_path(read_file)):
            genome_file = headers_of_genome_files[entry['target_description']]
            file_handles["%s-%s" % (read_file, genome_file)].write(
                segemehl_builder.entry_to_line(entry))
        for output_file in file_handles.values():
            output_file.close()

    def _headers_of_genome_files(self):
        """Extract the FASTA headers of all genome files."""
        headers = {}
        for genome_file in self.genome_files:
            genome_fh = open(self._genome_file_path(genome_file))
            headers[genome_fh.readline()[:-1]] = genome_file
        return(headers)

    def trace_reads_after_mapping(self):
        """Trace the way of each read during the different steps.

        A file is generated that can be used for downstream
        statistics.
        """
        for read_file in self.read_files:
            self.read_ids = []
            self.read_ids_and_traces = {}
            self._get_read_ids_and_lengths(read_file)
            self._read_first_mapping_output(read_file)
            self._read_first_mapping_unmapped_reads(read_file)
            self._read_clipped_unmapped_reads(read_file)
            self._read_size_filtered_reads_passed(read_file)
            self._read_size_filtered_reads_failed(read_file)
            self._read_second_mapping_output(read_file)
            self._read_second_mapping_unmapped_reads(read_file)
            self._read_combined_mapping_a_filtered_passed(read_file)
            self._read_combined_mapping_a_filtered_failed(read_file)
            self._write_trace_file(read_file)

    def _write_trace_file(self, read_file):
        """Write the trace of each read to a file.

        Arguments:
        - `read_file,`: the read file that was used to generate the
                        mapping files.

        """
        trace_fh = open(self._trace_file_path(read_file), "w")
        trace_fh.write("#Read id\tRead length\tNumber of mappings first run\t"
                       "length after clipping\tPassed size filter\t"
                       "Number of mappings second run\t"
                       "Passed a-content filter\tMapping length\t"
                       "Final status\n")
        for read_id in self.read_ids:
            trace = self.read_ids_and_traces[read_id]
            trace.setdefault("no_of_mappings_first_run", "-")
            trace.setdefault("length_after_clipping", "-")
            trace.setdefault("passed_size_filtering", "-")
            trace.setdefault("no_of_mappings_second_run", "-")
            trace.setdefault("passed_a-content_filtering", "-")
            trace.setdefault("mapping_length", "-")
            result_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                read_id, trace["length"], trace["no_of_mappings_first_run"],
                trace["length_after_clipping"], trace["passed_size_filtering"],
                trace["no_of_mappings_second_run"], 
                trace["passed_a-content_filtering"], 
                trace["mapping_length"],
                self._final_mapping_status(trace))
            trace_fh.write(result_line)

    def _get_read_ids_and_lengths(self, read_file):
        """Get the ids and length of read sequences.

        Here all the initial reads are indexed and the length
        saved. The created dictionary is filled further in following
        steps.

        Arguments:
        - `read_file`: name of the read file
        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self._read_file_path(read_file)):
            self.read_ids.append(header)
            self.read_ids_and_traces[header] = {'length' : len(seq)}

    def _read_first_mapping_output(self, read_file):
        """Read the result of the first Segemehl mapping.
        
        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the mapping file.

        """
        segemehl_parser = SegemehlParser()
        for entry in segemehl_parser.entries(
            self._raw_read_mapping_output_path(read_file)):
            entry_id = entry["id"][1:] # remove ">"
            self.read_ids_and_traces[entry_id].setdefault(
                "no_of_mappings_first_run", 0)
            self.read_ids_and_traces[entry_id]["no_of_mappings_first_run"] += 1

    def _read_first_mapping_unmapped_reads(self, read_file):
        """Read the file of unmapped reads of the first run.

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the mapping file.

        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self._unmapped_raw_read_file_path(read_file)):
                self.read_ids_and_traces[header]["no_of_mappings_first_run"] = 0

    def _read_clipped_unmapped_reads(self, read_file):
        """Read the file of clipped unmapped reads.

        Stores the length of the reads.

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the mapping file.

        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self._unmapped_read_clipped_path(read_file)):
            self.read_ids_and_traces[header][
                "length_after_clipping"] = len(seq)

    def _read_size_filtered_reads_passed(self, read_file):
        """Read the file of reads that passed the size filter.

        Arguments:
        - `read_file`: the name of the orignal read file

        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self._unmapped_clipped_size_filtered_read_path(read_file)):
            if header == "": continue
            self.read_ids_and_traces[header][
                "passed_size_filtering"] = True

    def _read_size_filtered_reads_failed(self, read_file):
        """Read the file of reads that did not pass the size filter.

        Arguments:
        - `read_file`: the name of the orignal read file

        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self._unmapped_clipped_size_failed_read_path(read_file)):
            if header == "": continue
            self.read_ids_and_traces[header][
                "passed_size_filtering"] = False

    def _read_second_mapping_output(self, read_file):
        """Read the file of the second Segemehl mapping.

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the mapping file.
        """
        segemehl_parser = SegemehlParser()
        for entry in segemehl_parser.entries(
            self._clipped_reads_mapping_output_path(read_file)):
            entry_id = entry["id"][1:] # remove ">"
            self.read_ids_and_traces[entry_id].setdefault(
                "no_of_mappings_second_run", 0)
            self.read_ids_and_traces[entry_id]["no_of_mappings_second_run"] += 1

    def _read_second_mapping_unmapped_reads(self, read_file):
        """Read the fasta file unmappable read of the second run.

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the first mapping file.
        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self._unmapped_reads_second_mapping_path(read_file)):
            if header == "": continue
            self.read_ids_and_traces[header]["no_of_mappings_second_run"] = 0

    def _read_combined_mapping_a_filtered_passed(self, read_file):
        """Read file of mappings passing the A-contend filtering.

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the first mapping file.
        """
        segemehl_parser = SegemehlParser()
        for entry in segemehl_parser.entries(
            self._combined_mapping_file_a_filtered_path(read_file)):
            entry_id = entry["id"][1:] # remove ">"
            self.read_ids_and_traces[entry_id][
                "passed_a-content_filtering"] = True
            self.read_ids_and_traces[entry_id][
                "mapping_length"] = len(entry["sequence"])
    
    def _read_combined_mapping_a_filtered_failed(self, read_file):
        """Read file of mapping not passing the A-contend filtering.

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the first mapping file.
        """
        segemehl_parser = SegemehlParser()
        for entry in segemehl_parser.entries(
            self._combined_mapping_file_a_filter_failed_path(read_file)):
            entry_id = entry["id"][1:] # remove ">"
            self.read_ids_and_traces[entry_id][
                "passed_a-content_filtering"] = False

    def build_gr_files(self):
        """Generate GR files for all read/genome file combinations."""
        for read_file in self.read_files:
            for genome_file in self.genome_files:
                self._build_gr_file(read_file, genome_file)

    def _build_gr_file(self, read_file, genome_file):
        """Generate GR files

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the first mapping file.
        - `genome_file`: name of the target genome file.
        """
        call("%s %s/segemehl2gr.py -o %s %s" % (
                self.python_bin, self.bin_folder,
                self._gr_file_path(read_file, genome_file),
                self._combined_mapping_file_a_filtered_split_path(read_file, 
                                                                  genome_file)),
             shell=True)

    def build_read_normalized_gr_files(self):
        """Generate normalized GR files for all read/genome files"""
        for genome_file in self.genome_files:
            lowest_number_of_mappings = self._lowest_number_of_mappings(
                genome_file)
            for read_file in self.read_files:
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
        call("%s %s/segemehl2gr.py -r -m %s -o %s %s" % (
                self.python_bin, self.bin_folder, lowest_number_of_mappings,
                self._gr_read_normalized_file_path(read_file, genome_file),
                self._combined_mapping_file_a_filtered_split_path(
                    read_file, genome_file)), shell=True)

    def build_nucl_normalized_gr_files(self):
        """Generate normalized GR files for all read/genome files"""
        for genome_file in self.genome_files:
            lowest_number_of_mapped_nucleotides = (
                self._lowest_number_of_mapped_nucleotides(genome_file))
            for read_file in self.read_files:
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
        call("%s %s/segemehl2gr.py -n -m %s -o %s %s" % (
                self.python_bin, self.bin_folder, 
                lowest_number_of_mapped_nucleotides,
                self._gr_nucl_normalized_file_path(read_file, genome_file),
                self._combined_mapping_file_a_filtered_split_path(
                    read_file, genome_file)), shell=True)

    def _lowest_number_of_mappings(self, genome_file):
        """Return the lowest number of mappings found.

        Arguments:
        - `genome_file`: target genome file

        """
        lowest_number_of_mappings = min(
            [self._count_mapped_reads(read_file, genome_file) 
             for read_file in self.read_files])
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
             for read_file in self.read_files])
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
        segemehl_parser = SegemehlParser()
        seen_ids = {}
        for entry in segemehl_parser.entries(
            self._combined_mapping_file_a_filtered_split_path(
                read_file, genome_file)):
            seen_ids[entry['id']] = 1
        return(len(seen_ids))
    
    def _count_mapped_nucleotides(self, read_file, genome_file):
        """Count number of successfully mapped reads.

        Arguments:
        - `read_file`: orignal read file used to generate the mappings.
        - `genome_file`: targe genome file

        """
        segemehl_parser = SegemehlParser()
        seen_ids = []
        nucleotide_counting = 0
        for entry in segemehl_parser.entries(
            self._combined_mapping_file_a_filtered_split_path(
                read_file, genome_file)):
            if entry['id'] in seen_ids:
                continue
            nucleotide_counting += len(entry['sequence'])
            seen_ids.append(entry['id'])
        return(nucleotide_counting)

    def find_annotation_hits(self):
        """Search for overlaps of reads and annotations."""
        for read_file in self.read_files:
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

    def _final_mapping_status(self, trace):
        """Return the final mapping status of a read.

        Arguments:
        - `trace`: the trace of the a read.
        """
        if (trace["passed_a-content_filtering"] and 
            not trace["passed_a-content_filtering"] == "-"):
            if trace["no_of_mappings_first_run"] > 0:
                return("mapped_in_first_round")
            elif trace["no_of_mappings_second_run"] > 0:
                return("mapped_in_second_round")
        elif not trace["passed_a-content_filtering"]:
            if trace["no_of_mappings_first_run"] > 0:
                return("mapped_in_first_round-failed_a-content_filter")
            elif trace["no_of_mappings_second_run"] > 0:
                return("mapped_in_second_round-faild_a_content_filter")
        elif trace["passed_a-content_filtering"] == "-":
            if not trace["passed_size_filtering"]:
                return("failed_size_filter_after_clipping")
            if trace["no_of_mappings_second_run"] == 0:
                return("not_mappable_in_second_run")
        else:
            return("lost_somewhere")

    def create_tracing_summay(self):
        """
        """
        stati = [
            "mapped_in_first_round", 
            "mapped_in_first_round-failed_a-content_filter",
            "mapped_in_second_round", 
            "mapped_in_second_round-faild_a_content_filter",
            "failed_size_filter_after_clipping", "not_mappable_in_second_run",
            "lost_somewhere"]
        summary_fh = open(self.tracing_summary_file, "w")
        summary_fh.write(
            "#lib name\t" + 
            "\t".join(stati) +
            "\t" +
            "total number of reads\t" +
            "sum of mappable reads\t" + 
            "percentage mappable reads\t" + 
            "\n")
        for read_file in self.read_files:
            stati_and_countings = self._summarize_tracing_file(
                self._trace_file_path(read_file))
            countings = []
            for status in stati:
                stati_and_countings.setdefault(status, 0)
                countings.append(stati_and_countings[status])
            summary_fh.write(
                read_file + "\t" + 
                "\t".join([str(counting) for counting in countings]) + 
                "\t%s" % sum(countings) +
                "\t%s" % (stati_and_countings["mapped_in_first_round"] + 
                          stati_and_countings["mapped_in_second_round"]) +
                "\t%s" % ((stati_and_countings["mapped_in_first_round"] + 
                          stati_and_countings["mapped_in_second_round"])/
                          sum(countings)*100.0) +
                "\n")
        summary_fh.close()

    def _summarize_tracing_file(self, tracing_file):
        """
        """
        stati_and_countings = {}
        for line in open(tracing_file):
            if line[0] in ["#", "\n"]:
                continue
            final_status = line[:-1].split("\t")[8]
            stati_and_countings.setdefault(final_status, 0)
            stati_and_countings[final_status] += 1
        return(stati_and_countings)

    def create_lib_genome_summary(self):
        """Create a file with lib/genome based read coutings."""
        coutings = []
        segemehl_parser = SegemehlParser()
        for read_file in self.read_files:
            counting_row = []
            for genome_file in self.genome_files:
                reads = {}
                for entry in segemehl_parser.entries(
                    self._combined_mapping_file_a_filtered_split_path(
                        read_file, genome_file)):
                    reads[entry["id"]] = 1
                counting_row.append(len(reads))
            coutings.append(counting_row)
        summary_fh = open(self.lib_genome_read_mapping_summary, "w")
        summary_fh.write("\t%s\n" % ("\t".join(self.genome_files)))
        for read_file, row in zip(self.read_files, coutings):
            summary_fh.write("%s\t%s\n" % (read_file, "\t".join(
                        [str(value) for value in row])))

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


        
# OBSOLETE

    # def _create_root_folder(self, project_name):
    #     """Create the root folder of a new project with the given name.
        
    #     Arguments:
    #     - `project_name`: Name of the project root folder

    #     """
    #     if not os.path.exists(project_name):
    #         os.mkdir(project_name)
    #     else:
    #         sys.stderr.write("Cannot create folder \"%s\"! File/folder with "
    #                          "the same name exists already.\n" % project_name)
    #         sys.exit(2)

    # def _create_subfolders(self, project_name):
    #     """Create required subfolders in the given folder.
        
    #     Arguments:
    #     - `project_name`: Name of the project root folder

    #     """
    #     for folder in [
    #         self.input_folder, self.output_folder, self.rna_seq_folder,
    #         self.annotation_folder, self.read_mappings_first_run_folder,
    #         self.read_mappings_second_run_folder, self.gr_folder,
    #         self.gr_folder_read_normalized, self.gr_folder_nucl_normalized,
    #         self.read_mapping_index_folder, self.genome_folder,
    #         self.umapped_reads_of_first_mapping_folder,
    #         self.umapped_reads_of_second_mapping_folder,
    #         self.combined_mappings_folder, self.combined_mapping_split_folder,
    #         self.annotation_hit_folder, self.annotation_hit_overview_folder,
    #         self.annotation_hit_overview_read_normalized_folder,
    #         self.annotation_hit_overview_nucl_normalized_folder,
    #         self.read_tracing_folder, self.input_file_stats_folder, 
    #         self.report_folder]:
    #         folder_in_root_folder = "%s/%s" % (project_name, folder)
    #         if not os.path.exists(folder_in_root_folder):
    #             os.mkdir(folder_in_root_folder)
    # def _create_config_file(self, project_name):
    #     """Create a JSON config file.
        
    #     Arguments:
    #     - `project_name`: Name of the project root folder
    #     """
    #     config_fh = open("%s/%s" % (project_name, self.config_file), "w")
    #     config_fh.write(json.dumps({"annotation_and_genomes_files" : {}}))
    #     config_fh.close()


    # def _create_read_file_stats(self):
    #     """Create a stat file for the input read files."""
    #     stat_fh = open(self.read_file_stats, "w")
    #     for read_file in self.read_files:
    #         stat_fh.write("%s:\n" % (read_file))
    #         stat_fh.write("* SHA256: %s\n" % (
    #                 self._sha256_of_file(self._read_file_path(read_file))))
    #         stat_fh.write("* Number of lines: %s\n" % (
    #                 self._number_of_lines_in_file(
    #                     self._read_file_path(read_file))))
    #         stat_fh.write("* Size: %s\n" % (
    #                 self._file_size(self._read_file_path(read_file))))
    #         stat_fh.write("* Number of Fasta entries: %s\n" % (
    #                 self._number_of_fasta_entries(
    #                     self._read_file_path(read_file))))
    #         stat_fh.write("\n")
    #     stat_fh.close()

    # def _create_genome_file_stats(self):
    #     """Create a stat file for the input genome files."""
    #     stat_fh = open(self.genome_file_stats, "w")
    #     for genome_file in self.genome_files:
    #         stat_fh.write("%s:\n" % (genome_file))
    #         stat_fh.write("* SHA256: %s\n" % (self._sha256_of_file(
    #                     self._genome_file_path(genome_file))))
    #         stat_fh.write("* Number of lines: %s\n" % (
    #                 self._number_of_lines_in_file(
    #                     self._genome_file_path(genome_file))))
    #         stat_fh.write("* Size: %s\n" % (self._file_size(
    #                                 self._genome_file_path(genome_file))))
    #         stat_fh.write("\n")
    #     stat_fh.close()

    # def _sha256_of_file(self, file_path):
    #     """Calculate the SHA256 hash sum of a given file

    #     Arguments:
    #     - `file_path`: path of the file to process

    #     """
    #     # Todo: Fix handling of large files and then use this again
    #     # instead of calling the shell command
    #     #return(hashlib.sha256(open(file_path).read().encode()).hexdigest())
    #     return(Popen("sha256sum %s" % file_path , shell=True, stdout=PIPE
    #                  ).communicate()[0].split()[0].decode("utf-8"))

    # def _file_size(self, file_path):
    #     """Return the size of a given file.

    #     Arguments:
    #     - `file_path`: path of the file to process

    #     """
    #     return(os.path.getsize(file_path))

    # def _number_of_lines_in_file(self, file_path):
    #     """Return the number of lines of a given file.

    #     Arguments:
    #     - `file_path`: path of the file to process

    #     """
    #     return(len(open(file_path).readlines()))

    # def _number_of_fasta_entries(self, file_path):
    #     """Return the number of fasta entries.

    #     Arguments:
    #     - `file_path`: path of the file to process
        
    #     """
    #     fasta_parser = FastaParser()
    #     counter = 0
    #     for head, seq in fasta_parser.parse_fasta_file(file_path):
    #         counter += 1
    #     return(counter)
