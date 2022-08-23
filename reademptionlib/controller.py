import concurrent.futures
import os
import sys
import json
import pysam
from reademptionlib.bamsorter import BamSorter
from reademptionlib.coveragecreator import CoverageCreator
from reademptionlib.crossalignfilter import CrossAlignFilter
from reademptionlib.deseq import DESeqRunner
from reademptionlib.fasta import FastaParser
from reademptionlib.genewisequanti import GeneWiseQuantification
from reademptionlib.genewisequanti import GeneWiseOverview
from reademptionlib.pathcreator import PathCreator
from reademptionlib.projectcreator import ProjectCreator
from reademptionlib.rawstatdata import RawStatDataWriter, RawStatDataReader
from reademptionlib.readaligner import ReadAligner
from reademptionlib.readalignerstats import ReadAlignerStats
from reademptionlib.readalignerstatstable import ReadAlignerStatsTable
from reademptionlib.readprocessor import ReadProcessor
from reademptionlib.vizalign import AlignViz
from reademptionlib.vizdeseq import DESeqViz
from reademptionlib.vizgenequanti import GeneQuantiViz
from reademptionlib.fragmentbuilder import FragmentBuilder
from datetime import datetime

class Controller(object):

    """Manage the actions of the subcommands.

    The Controller takes care of providing the arguments like path
    names and the parallel processing of tasks.

    """

    def __init__(self, args):
        """Create an instance."""
        self._args = args
        # the species argument containing the information about species-folder
        # naming and species display-names is only provided by the user when
        # creating a new project via the argument '--species'.
        # When creating a new project the information is written into the
        # config file.
        # For all other subcommands the information is not provided via an
        # argument but will be retrieved by path_creator from the
        # config file
        if "species" not in self._args:
            self._args.species = None
        self._pathcreator = PathCreator(args.project_path, self._args.species)
        self._species_folder_prefixes_and_display_names = (
            self._pathcreator.species_folder_prefixes_and_display_names
        )
        self._species_folder_prefixes = (
            self._pathcreator.species_folder_prefixes
        )
        # Some words are used for creating specific folders or files.
        # These words can not be used as species folder or display names.
        # If one word of the forbidden list below is used READemption exits and
        # will not create any folders
        self.forbidden_species_folder_prefixes = ["all", "read_lengths"]
        self._species_display_names = self._pathcreator.species_display_names
        if (
            not set(self.forbidden_species_folder_prefixes).isdisjoint(
                set(self._species_folder_prefixes)
            )
        ) or (
            not set(self.forbidden_species_folder_prefixes).isdisjoint(
                set(self._species_display_names)
            )
        ):
            error_message = (
                "The following words are not allowed as a SUB_FOLDER_NAME or "
                "DISPLAY_NAME: "
                f"{self.forbidden_species_folder_prefixes} \n"
                "Please rerun the command with SUB_FOLDER_NAMES and "
                "DISPLAY_NAMES that are not in the list. \n"
            )
            self._write_err_msg_and_quit(error_message)

        self._read_files = None
        self._ref_seq_files = None

    def create_project(self, version):
        """Create a new project."""
        sys.stdout.write(
            "   ___  _______   ___                 __  _\n"
            "  / _ \\/ __/ _ | / _ \\___ __ _  ___  / /_(_)__  ___\n"
            " / , _/ _// __ |/ // / -_)  ' \\/ _ \\/ __/ / _ \\/ _ \\\n"
            "/_/|_/___/_/ |_/____/\\__/_/_/_/ .__/\\__/_/\\___/_//_/\n"
            "                             / /\n"
            "====================================================\n"
            "========================================\n"
            "=======================\n"
            "==============\n\n"
            "[https://reademption.readthedocs.io/en/latest/]\n\n"
        )
        project_creator = ProjectCreator()
        project_creator.create_root_folder(self._args.project_path)
        project_creator.create_config_file(
            self._pathcreator.config_file,
            self._species_folder_prefixes_and_display_names,
        )
        project_creator.create_subfolders(
            self._pathcreator.required_new_project_folders()
        )
        project_creator.create_version_file(
            self._pathcreator.version_path, version
        )
        sys.stdout.write(
            'Created folder "%s" and required subfolders.\n'
            % (self._args.project_path)
        )
        ref_seq_folders = ", ".join(
            (
                f'"{folder}"'
                for folder in self._pathcreator.ref_seq_folders_by_species.values()
            )
        )
        sys.stdout.write(
            f'Please copy read files into folder "{self._pathcreator.read_fasta_folder}" and '
            f"reference sequences files into folder/s {ref_seq_folders}.\n"
        )

    def _get_references_by_species(self) -> dict:
        """
        Reads the reference sequence fasta files in the species folders and
        creates a dictionary containing the species as keys and a list
        containing the reference ids for each species as values:
        Example:
        {'human': ['GL000008.2', 'chr1', 'chr2', 'chr3', 'chr19', 'chr20',
                   'chr21', 'chr22', 'chr23', 'chrX', 'KQ031388.1'],
        'staphylococcus': ['NC_007795.1'],
        'influenza': ['NC_007373.1', 'NC_007372.1']}
        :return: the reference ids for each species
        """
        fasta_parser = FastaParser()
        references_by_species = {}
        for (
            species,
            reference_files,
        ) in self._pathcreator.ref_seq_paths_by_species.items():
            references_by_species[species] = []
            for reference_file in reference_files:
                with open(reference_file, "r") as fasta_fh:
                    for header, sequence in fasta_parser.entries(fasta_fh):
                        header_id = fasta_parser.header_id(header)
                        references_by_species[species].append(header_id)
        return references_by_species

    def align_reads(self):
        """Perform the alignment of the reads."""
        self._test_folder_existance(
            self._pathcreator.required_read_alignment_folders()
        )
        assert self._args.paired_end in [True, False]
        self._pathcreator.set_ref_seq_paths_by_species()
        self._ref_seq_files = self._pathcreator.get_ref_seq_files()
        self._pathcreator.set_ref_seq_path_list()
        self._test_align_file_existance()
        if not self._args.paired_end:
            # Single end reads
            self._read_files = self._pathcreator.get_read_files()
            self._lib_names = self._pathcreator.get_lib_names_single_end()
            self._pathcreator.set_read_files_dep_file_lists_single_end(
                self._read_files, self._lib_names
            )
            self._prepare_reads_single_end()
            print(f"controller align_single_end_reads start {datetime.now()}")
            self._align_single_end_reads()
            print(f"controller align_single_end_reads stop {datetime.now()}")
        else:
            # Paired end reads
            self._read_file_pairs = self._pathcreator.get_read_file_pairs()
            self._lib_names = self._pathcreator.get_lib_names_paired_end()
            self._pathcreator.set_read_files_dep_file_lists_paired_end(
                self._read_file_pairs, self._lib_names
            )
            self._prepare_reads_paired_end()
            print(f"controller align_paired_end_reads start {datetime.now()}")
            self._align_paired_end_reads()
            print(f"controller align_paired_end_reads stop {datetime.now()}")
        print(f"controller generate_read_alignment_stats start {datetime.now()}")
        self._generate_read_alignment_stats(
            self._lib_names,
            self._pathcreator.read_alignment_bam_paths,
            self._pathcreator.unaligned_reads_paths,
            self._pathcreator.read_alignments_stats_path,
            self._args.paired_end
        )
        print(f"controller generate_read_alignment_stats stop {datetime.now()}")
        if self._args.crossalign_cleaning:
            self._remove_crossaligned_reads()

        if self._args.paired_end:
            # Build a bam file containing fragments merged from read
            # pairs
            if not self._args.no_fragment_building:
                fragments = True
                # sort the bam files by name and sam tag hit index to
                # accelerate fragment building
                print(f"controller sort bams by name and index start {datetime.now()}")
                self._sort_bams_by_name_and_index()
                print(f"controller sort bams by name and index end {datetime.now()}")
                # build the fragments bam file
                print(f"controller build_fragments start {datetime.now()}")
                self._build_fragments()
                print(f"controller build_fragments stop {datetime.now()}")
                # generate fragment alignment stats
                print(f"controller generate_fragment_alignmnet_stats start {datetime.now()}")
                self._generate_read_alignment_stats(
                    self._lib_names,
                    self._pathcreator.aligned_fragments_bam_paths,
                    self._pathcreator.unaligned_reads_paths,
                    self._pathcreator.fragment_alignments_stats_path,
                    self._args.paired_end,
                    fragments
                )
                print(f"controller generate_fragment_alignmnet_stats stop {datetime.now()}")
                # write fragment stats table
                print(f"controller write_alignment_stats_table fragments start {datetime.now()}")
                self._write_alignment_stat_table(self._pathcreator.fragment_alignments_stats_path,
                                                 self._pathcreator.fragment_alignment_stats_table_path,
                                                 self._pathcreator.fragment_alignment_stats_table_transposed_path,
                                                 fragments)
                print(f"controller write_alignment_stats_table fragments stop {datetime.now()}")
        print(f"controller write_alignment_stats_table reads start {datetime.now()}")
        self._write_alignment_stat_table(self._pathcreator.read_alignments_stats_path,
                                         self._pathcreator.read_alignment_stats_table_path,
                                         self._pathcreator.read_alignment_stats_table_transposed_path)
        print(f"controller write_alignment_stats_table reads stop {datetime.now()}")

    def _sort_bams_by_name_and_index(self):
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=self._args.processes
        ) as executor:
            for (
                    read_alignment_path,
                    read_alignment_sorted_path,
            ) in zip(
                self._pathcreator.read_alignment_bam_paths,
                self._pathcreator.read_alignment_bam_sorted_paths,
            ):
                # Sort fragments
                bam_sorter = BamSorter("HI", sort_by="name")
                jobs.append(
                    executor.submit(bam_sorter.sort_bam,
                                    read_alignment_path,
                                    read_alignment_sorted_path)
                )
        # Evaluate thread outcome
        self._check_job_completeness(jobs)

    def _build_fragments(self):
        # Build a bam file containing fragments merged from read
        # pairs
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=self._args.processes
        ) as executor:
            for (
                    read_alignment_sorted_path,
                    fragment_alignment_path,
            ) in zip(
                self._pathcreator.read_alignment_bam_sorted_paths,
                self._pathcreator.aligned_fragments_bam_paths,
            ):
                # Perform the building for fragments from reads
                fragment_builder = FragmentBuilder(self._args.max_fragment_length)
                jobs.append(
                    executor.submit(fragment_builder.build_bam_file_with_fragments,
                                    read_alignment_sorted_path,
                                    fragment_alignment_path)
                )
            # Evaluate thread outcome
        self._check_job_completeness(jobs)


    def _remove_crossaligned_reads(self):
        # self._string_to_species_and_sequence_ids()
        references_by_species = self._get_references_by_species()
        # Determine cross-mapped reads
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes
        ) as executor:
            for (
                bam_path,
                bam_with_crossmappings_path,
                bam_cleaned_tmp_path,
                crossmapped_reads_path,
            ) in zip(
                self._pathcreator.read_alignment_bam_paths,
                self._pathcreator.read_alignment_bam_with_crossmappings_paths,
                self._pathcreator.read_alignment_bam_cross_cleaned_tmp_paths,
                self._pathcreator.crossmapped_reads_paths,
            ):
                # Perform the removal or cross aligned reads
                cross_align_filter = CrossAlignFilter(
                    bam_path,
                    bam_cleaned_tmp_path,
                    crossmapped_reads_path,
                    references_by_species,
                )
                jobs.append(
                    executor.submit(cross_align_filter.filter_crossmapped_reads)
                )
        # Evaluate thread outcome
        self._check_job_completeness(jobs)

        for (
            bam_path,
            bam_with_crossmappings_path,
            bam_cleaned_tmp_path,
            crossmapped_reads_path,
        ) in zip(
            self._pathcreator.read_alignment_bam_paths,
            self._pathcreator.read_alignment_bam_with_crossmappings_paths,
            self._pathcreator.read_alignment_bam_cross_cleaned_tmp_paths,
            self._pathcreator.crossmapped_reads_paths,
        ):
            # Rename the original mapping file that potentially
            # contains cross aligned reads
            os.rename(bam_path, bam_with_crossmappings_path)
            os.rename(bam_path + ".bai", bam_with_crossmappings_path + ".bai")
            # Move the cross aligned filtered file to the final mapping
            # path
            os.rename(bam_cleaned_tmp_path, bam_path)
            os.rename(bam_cleaned_tmp_path + ".bai", bam_path + ".bai")

    def _test_align_file_existance(self):
        """Test if the input file for the the align subcommand exist."""
        if len(self._pathcreator.get_read_files()) == 0:
            self._write_err_msg_and_quit("Error! No read libraries given!\n")
        if len(self._ref_seq_files) == 0:
            self._write_err_msg_and_quit(
                "Error! No reference sequence files given!\n"
            )

    def _test_folder_existance(self, task_specific_folders):
        """Test the existance of required folders."""
        for folder in (
            self._pathcreator.required_base_folders() + task_specific_folders
        ):
            if not os.path.exists(folder):
                self._write_err_msg_and_quit(
                    "Error! Folder '%s' does not exist! Is the given project "
                    "folder name correct?\n" % folder
                )

    def _file_needs_to_be_created(self, file_path, quiet=False):
        """Test if a file exists of need to be created."""
        if not self._args.check_for_existing_files:
            return True
        if os.path.exists(file_path):
            if not quiet:
                sys.stderr.write(
                    "File %s exists. Skipping its generation.\n" % file_path
                )
            return False
        return True

    def _prepare_reads_single_end(self):
        """Manage the prepartion of reads before the actual mappings."""
        read_files_and_jobs = {}
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes
        ) as executor:
            for lib_name, read_path, processed_read_path in zip(
                self._lib_names,
                self._pathcreator.read_paths,
                self._pathcreator.processed_read_paths,
            ):
                if not self._file_needs_to_be_created(processed_read_path):
                    continue
                read_processor = ReadProcessor(
                    poly_a_clipping=self._args.poly_a_clipping,
                    min_read_length=self._args.min_read_length,
                    fastq=self._args.fastq,
                    min_phred_score=self._args.min_phred_score,
                    adapter=self._args.adapter,
                    reverse_complement=self._args.reverse_complement,
                )
                read_files_and_jobs[lib_name] = executor.submit(
                    read_processor.process_single_end,
                    read_path,
                    processed_read_path,
                )
        self._evaluet_job_and_generate_stat_file(read_files_and_jobs)

    def _prepare_reads_paired_end(self):
        read_files_and_jobs = {}
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes
        ) as executor:
            for lib_name, read_path_pair, processed_read_path_pair in zip(
                self._lib_names,
                self._pathcreator.read_path_pairs,
                self._pathcreator.processed_read_path_pairs,
            ):
                for processed_read_path in processed_read_path_pair:
                    if not self._file_needs_to_be_created(processed_read_path):
                        continue
                    read_processor = ReadProcessor(
                        poly_a_clipping=False,
                        min_read_length=self._args.min_read_length,
                        fastq=self._args.fastq,
                        min_phred_score=self._args.min_phred_score,
                        adapter=self._args.adapter,
                        reverse_complement=self._args.reverse_complement,
                    )
                    read_files_and_jobs[lib_name] = executor.submit(
                        read_processor.process_paired_end,
                        read_path_pair,
                        processed_read_path_pair,
                    )
        self._evaluet_job_and_generate_stat_file(read_files_and_jobs)

    def _evaluet_job_and_generate_stat_file(self, read_files_and_jobs):
        raw_stat_data_writer = RawStatDataWriter(pretty=True)
        # Evaluate thread outcome
        self._check_job_completeness(read_files_and_jobs.values())
        if not self._file_needs_to_be_created(
            self._pathcreator.read_processing_stats_path
        ):
            return
        # Create a dict of the read file names and the processing
        # counting results
        read_files_and_stats = dict(
            [
                (lib_name, job.result())
                for lib_name, job in read_files_and_jobs.items()
            ]
        )
        raw_stat_data_writer.write(
            read_files_and_stats, self._pathcreator.read_processing_stats_path
        )

    def _align_single_end_reads(self):
        """Manage the actual alignment of single end reads."""
        read_aligner = ReadAligner(self._args.segemehl_bin, self._args.progress)
        if self._file_needs_to_be_created(self._pathcreator.index_path):
            read_aligner.build_index(
                self._pathcreator.ref_seq_path_list,
                self._pathcreator.index_path,
            )
        for read_path, output_path, nomatch_path in zip(
            self._pathcreator.processed_read_paths,
            self._pathcreator.read_alignment_bam_paths,
            self._pathcreator.unaligned_reads_paths,
        ):
            if not self._file_needs_to_be_created(output_path):
                continue

            read_aligner.run_alignment(
                read_path,
                self._pathcreator.index_path,
                self._pathcreator.ref_seq_path_list,
                output_path,
                nomatch_path,
                int(self._args.processes),
                int(self._args.segemehl_accuracy),
                float(self._args.segemehl_evalue),
                self._args.split,
                paired_end=False,
            )

    def _align_paired_end_reads(self):
        """Manage the actual alignemnt of paired end reads."""
        read_aligner = ReadAligner(self._args.segemehl_bin, self._args.progress)
        if self._file_needs_to_be_created(self._pathcreator.index_path):
            read_aligner.build_index(
                self._pathcreator.ref_seq_path_list,
                self._pathcreator.index_path,
            )
        for read_path_pair, output_path, nomatch_path in zip(
            self._pathcreator.processed_read_path_pairs,
            self._pathcreator.read_alignment_bam_paths,
            self._pathcreator.unaligned_reads_paths,
        ):
            if not self._file_needs_to_be_created(output_path):
                continue
            read_aligner.run_alignment(
                read_path_pair,
                self._pathcreator.index_path,
                self._pathcreator.ref_seq_path_list,
                output_path,
                nomatch_path,
                int(self._args.processes),
                int(self._args.segemehl_accuracy),
                float(self._args.segemehl_evalue),
                self._args.split,
                paired_end=True,
            )

    def _generate_read_alignment_stats(
        self,
        lib_names,
        result_bam_paths,
        unaligned_reads_paths,
        output_stats_path,
        paired_end=False,
        fragments=False
    ):
        """Manage the generation of alingment statistics."""
        raw_stat_data_writer = RawStatDataWriter(pretty=True)
        references_by_species = self._get_references_by_species()
        read_files_and_jobs = {}
        if not self._file_needs_to_be_created(output_stats_path):
            return
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes
        ) as executor:
            for (
                lib_name,
                read_alignment_bam_path,
                unaligned_reads_path,
            ) in zip(lib_names, result_bam_paths, unaligned_reads_paths):
                read_aligner_stats = ReadAlignerStats(references_by_species, paired_end, fragments)
                read_files_and_jobs[lib_name] = executor.submit(
                    read_aligner_stats.count,
                    read_alignment_bam_path,
                    unaligned_reads_path,
                )
        # Evaluate thread outcome
        self._check_job_completeness(read_files_and_jobs.values())
        read_files_and_stats = dict(
            [
                (lib_name, job.result())
                for lib_name, job in read_files_and_jobs.items()
            ]
        )
        raw_stat_data_writer.write(read_files_and_stats, output_stats_path)

    def _write_alignment_stat_table(self, alignments_stats_path, alignment_stats_table_path, alignment_stats_table_transposed_path, fragments=False):
        """Manage the creation of the mapping statistic output table."""
        raw_stat_data_reader = RawStatDataReader()
        read_processing_stats = raw_stat_data_reader.read(
            self._pathcreator.read_processing_stats_path
        )
        final_alignment_stats = raw_stat_data_reader.read(
            alignments_stats_path
        )
        references_by_species = self._get_references_by_species()
        read_aligner_stats_table = ReadAlignerStatsTable(
            read_processing_stats,
            final_alignment_stats,
            self._lib_names,
            alignment_stats_table_path,
            alignment_stats_table_transposed_path,
            self._args.paired_end,
            self._species_folder_prefixes_and_display_names,
            references_by_species,
            fragments
        )
        read_aligner_stats_table.write()

    def create_coverage_files(self):
        """Create coverage files based on the read alignments.

        The coverages are calculated per replicon and the results are
        written to the output file. This might be slower but if all
        coverages are determined at once the data structure will become
        too large when working with large reference sequences.

        """

        # Select if normalisation is based on fragment numbers
        if self._args.paired_end and not self._args.no_norm_by_fragments:
            norm_by_fragments = True
        else:
            norm_by_fragments = False
        if norm_by_fragments:
            reads_or_fragments = "fragments"
            alignment_stats_path = self._pathcreator.fragment_alignments_stats_path
        else:
            reads_or_fragments = "reads"
            alignment_stats_path = self._pathcreator.read_alignments_stats_path

        # Get alignment stats
        raw_stat_data_reader = RawStatDataReader()
        alignment_stats = [
            raw_stat_data_reader.read(
                alignment_stats_path
            )
        ]
        # Lib names was paired end
        lib_names = list(alignment_stats[0].keys())
        was_paired_end_alignment = self._was_paired_end_alignment(lib_names)

        # Quit if the wrong parameters have been chosen for the subcommand
        if (was_paired_end_alignment and not self._args.paired_end):
            self._write_err_msg_and_quit("The alignemnt seems to be based on paired end reads. "
                                         "Please also set "
                                         "the option '-P' or '--paired_end'.\n")

        if (self._args.no_fragments and not self._args.paired_end):
            self._write_err_msg_and_quit("The option '-nf' or "
                                         "'--no_fragments' is only valid "
                                         "for paired end reads. If you have "
                                         "paired end reads, please also set "
                                         "the option '-P' or '--paired_end'.\n")

        if (self._args.no_norm_by_fragments and not self._args.paired_end):
            self._write_err_msg_and_quit("The option '-nnf' or "
                                         "'--no_norm_by_fragments' is only valid "
                                         "for paired end reads. If you have "
                                         "paired end reads, please also set "
                                         "the option '-P' or '--paired_end'.\n")



        # Set read files and lib names
        if not was_paired_end_alignment:
            self._pathcreator.set_read_files_dep_file_lists_single_end(
                self._pathcreator.get_read_files(), lib_names
            )
        else:
            self._pathcreator.set_read_files_dep_file_lists_paired_end(
                self._pathcreator.get_read_files(), lib_names
            )
        # If fragments should be used and they were not created during alignment,
        # they will be created now
        if (not self._args.no_fragments and self._args.paired_end):
            bam_files_exist = []
            for bam_fragment_path in self._pathcreator.aligned_fragments_bam_paths:
                bam_files_exist.append(os.path.exists(bam_fragment_path))
            # If any of the bam files containing fragments is missing, create all
            # of them
            if not all(bam_files_exist):
                self._build_fragments()

        # Set alignment paths to fragments or single reads
        if (not self._args.no_fragments and self._args.paired_end):
            alignment_paths = self._pathcreator.aligned_fragments_bam_paths
        else:
            alignment_paths = self._pathcreator.read_alignment_bam_paths
        # determine species cross mapped reads
        self._pathcreator.set_ref_seq_paths_by_species()
        if not self._args.count_cross_aligned_reads:
            self._crossmapped_reads_by_lib = {}
            for lib_name, read_alignment_path in zip(
                    lib_names, self._pathcreator.read_alignment_bam_paths
            ):
                # retrieve the cross mapped reads from the single read files
                # to also get reads where two mates map to different
                # species. This would not be possible with the built fragments
                self._crossmapped_reads_by_lib[
                    lib_name
                ] = self.determine_crossmapped_reads(read_alignment_path)

        if not self._args.non_strand_specific:
            strands = ["forward", "reverse"]
        else:
            strands = ["forward_and_reverse"]



        self.read_files_aligned_read_freq_and_min_reads_aligned_by_species = {}
        for sp in self._species_folder_prefixes_and_display_names.keys():
            # Retrieve the either the no. of uniquely aligned reads or
            # the number of species exclusive aligned reads ("all aligned" - "cross aligned") (Default behaviour) or
            # number of all aligned reads for each library of the given species
            read_files_aligned_read_freq = {}
            for read_file, attributes in alignment_stats[0].items():
                # If option normalize by uniquely is chosen, only the sum of uniquely aligned reads is used for normalisation
                # this excludes species cross mapped reads, split aligned reads and multiple aligned reads
                if self._args.normalize_by_uniquely:
                    read_files_aligned_read_freq[read_file] = attributes[
                        "species_stats"
                    ][sp][f"no_of_uniquely_aligned_{reads_or_fragments}"]
                elif self._args.normalize_cross_aligned_reads_included:
                    read_files_aligned_read_freq[read_file] = attributes[
                        "species_stats"
                    ][sp][f"no_of_aligned_{reads_or_fragments}"]
                # Default: Number of aligned reads without the cross aligned reads are used for normalization
                else:
                    read_files_aligned_read_freq[read_file] = (
                        attributes["species_stats"][sp][f"no_of_aligned_{reads_or_fragments}"]
                        - attributes["species_stats"][sp][
                            f"no_of_cross_aligned_{reads_or_fragments}"
                        ]
                    )
            self.read_files_aligned_read_freq_and_min_reads_aligned_by_species[
                sp
            ] = {}
            self.read_files_aligned_read_freq_and_min_reads_aligned_by_species[
                sp
            ]["read_files_aligned_read_freq"] = read_files_aligned_read_freq
            # Retrieve the min no. of aligned reads
            # of all libraries for the given species
            min_no_of_aligned_reads = float(
                min(read_files_aligned_read_freq.values())
            )
            self.read_files_aligned_read_freq_and_min_reads_aligned_by_species[
                sp
            ]["min_no_of_aligned_reads"] = min_no_of_aligned_reads
        self._pathcreator.set_coverage_folder_and_file_names(
            strands,
            lib_names,
            self.read_files_aligned_read_freq_and_min_reads_aligned_by_species,
        )

        project_creator = ProjectCreator()
        project_creator.create_subfolders(
            self._pathcreator.required_coverage_folders()
        )
        self._test_folder_existance(
            self._pathcreator.required_coverage_folders()
        )



        # get references by species
        references_by_species = self._get_references_by_species()

        # Run the coverage file creation species-wise
        for (
            sp
        ) in (
            self.read_files_aligned_read_freq_and_min_reads_aligned_by_species.keys()
        ):
            # Run the generation of coverage in parallel

            jobs = []
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=self._args.processes
            ) as executor:
                for lib_name, bam_path in zip(
                    lib_names, alignment_paths
                ):
                    if not self._args.count_cross_aligned_reads:
                        cross_mapped_reads = self._crossmapped_reads_by_lib[
                            lib_name
                        ]
                    else:
                        cross_mapped_reads = None
                    coverage_creator = CoverageCreator(
                        self._args,
                        strands,
                        self._pathcreator.coverage_files_by_species[sp][
                            lib_name
                        ],
                        references_by_species[sp],
                        self._args.count_cross_aligned_reads,
                        cross_mapped_reads,
                    )
                    jobs.append(
                        executor.submit(
                            coverage_creator.create_coverage_files_for_lib,
                            lib_name,
                            bam_path,
                            self.read_files_aligned_read_freq_and_min_reads_aligned_by_species[
                                sp
                            ][
                                "read_files_aligned_read_freq"
                            ][
                                lib_name
                            ],
                            self.read_files_aligned_read_freq_and_min_reads_aligned_by_species[
                                sp
                            ][
                                "min_no_of_aligned_reads"
                            ],
                        )
                    )
            # Evaluate thread outcome
            self._check_job_completeness(jobs)

    def _check_job_completeness(self, jobs):
        """Check the completness of each job in a list"""
        for job in concurrent.futures.as_completed(jobs):
            if job.exception():
                raise (job.exception())

    def determine_crossmapped_reads(self, read_alignment_path):
        """Find reads that are mapped to sequences of different species.

        The comparison is performed replicon wise in order to reduce
        the memory footprint.
        """
        references_by_species = self._get_references_by_species()
        crossmapped_reads = set()
        done_replicon_comparison = []
        with pysam.AlignmentFile(read_alignment_path) as bam:
            for org, replicon_ids in references_by_species.items():
                for replicon_id in replicon_ids:
                    self._read_ids = set()
                    # First, collect the ids of the aligned reads of
                    # this replicon
                    for alignment in bam.fetch(reference=replicon_id):
                        self._read_ids.add(alignment.qname)
                    # Then compare them to the alignments of each
                    # replicon of the other organism(s)
                    for (
                        comp_org,
                        comp_replicon_ids,
                    ) in references_by_species.items():
                        # Only compare replicons of different species
                        if org == comp_org:
                            continue
                        for comp_replicon_id in comp_replicon_ids:
                            comparison = sorted([replicon_id, comp_replicon_id])
                            # Check if comparison of the two replicons
                            # has been done already
                            if comparison in done_replicon_comparison:
                                continue
                            done_replicon_comparison.append(comparison)
                            # Compare all read ids of the comparison
                            # replicon to the query replicon read ids
                            for alignment in bam.fetch(
                                reference=comp_replicon_id
                            ):
                                if alignment.qname in self._read_ids:
                                    crossmapped_reads.add(alignment.qname)
        no_of_crossmapped_reads = len(crossmapped_reads)
        return crossmapped_reads

    def quantify_gene_wise(self):
        """Manage the counting of aligned reads per gene."""
        # Get alignment stats
        raw_stat_data_reader = RawStatDataReader()
        alignment_stats = [
            raw_stat_data_reader.read(
                self._pathcreator.read_alignments_stats_path
            )
        ]

        # Lib names and was paired end
        lib_names = sorted(list(alignment_stats[0].keys()))
        was_paired_end_alignment = self._was_paired_end_alignment(lib_names)

        # Quit if the wrong parameters have been chosen for the subcommand
        if (was_paired_end_alignment and not self._args.paired_end):
            self._write_err_msg_and_quit("The alignemnt seems to be based on paired end reads. "
                                         "Please also set "
                                         "the option '-P' or '--paired_end'.\n")

        if (self._args.no_fragments and not self._args.paired_end):
            self._write_err_msg_and_quit("The option '-nf' or "
                                         "'--no_fragments' is only valid "
                                         "for paired end reads. If you have "
                                         "paired end reads, please also set "
                                         "the option '-P' or '--paired_end'.\n")

        if (self._args.no_norm_by_fragments and not self._args.paired_end):
            self._write_err_msg_and_quit("The option '-nnf' or "
                                         "'--no_norm_by_fragments' is only valid "
                                         "for paired end reads. If you have "
                                         "paired end reads, please also set "
                                         "the option '-P' or '--paired_end'.\n")
        # Select if normalisation is based on fragment numbers
        if self._args.paired_end and not self._args.no_norm_by_fragments:
            norm_by_fragments = True
        else:
            norm_by_fragments = False

        project_creator = ProjectCreator()
        project_creator.create_subfolders(
            self._pathcreator.required_gene_quanti_folders()
        )
        self._test_folder_existance(
            self._pathcreator.required_gene_quanti_folders()
        )
        norm_by_alignment_freq = True
        norm_by_overlap_freq = True
        self._pathcreator.set_ref_seq_paths_by_species()
        references_by_species = self._get_references_by_species()
        if self._args.no_count_split_by_alignment_no:
            norm_by_alignment_freq = False
        if self._args.no_count_splitting_by_gene_no:
            norm_by_overlap_freq = False
        self._pathcreator.set_annotation_paths_by_species()

        # Set read files and lib names
        if not was_paired_end_alignment:
            self._pathcreator.set_read_files_dep_file_lists_single_end(
                self._pathcreator.get_read_files(), lib_names
            )
        else:
            self._pathcreator.set_read_files_dep_file_lists_paired_end(
                self._pathcreator.get_read_files(), lib_names
            )

        # If fragments should be used and they were not created during alignment,
        # they will be created now
        if (not self._args.no_fragments and self._args.paired_end):
            bam_files_exist = []
            for bam_fragment_path in self._pathcreator.aligned_fragments_bam_paths:
                bam_files_exist.append(os.path.exists(bam_fragment_path))
            # If any of the bam files containing fragments is missing, create all
            # of them
            if not all(bam_files_exist):
                self._build_fragments()

        # Set alignment paths to fragments or single reads
        self._pathcreator.set_annotation_files_by_species()
        if (not self._args.no_fragments and self._args.paired_end):
            alignment_paths = self._pathcreator.aligned_fragments_bam_paths
        else:
            alignment_paths = self._pathcreator.read_alignment_bam_paths

        # determine species cross mapped reads
        if not self._args.count_cross_aligned_reads:
            self._crossmapped_reads_by_lib = {}
            for lib_name, read_alignment_path in zip(
                lib_names, self._pathcreator.read_alignment_bam_paths
            ):
                # retrieve the cross mapped reads from the single read files
                # to also get reads where two mates map to different
                # species. This would not be possible with the built fragments
                self._crossmapped_reads_by_lib[
                    lib_name
                ] = self.determine_crossmapped_reads(read_alignment_path)
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._args.processes
        ) as executor:
            for (
                sp,
                annotation_files,
            ) in self._pathcreator.annotation_files_by_species.items():
                for lib_name, read_alignment_path in zip(
                    lib_names, alignment_paths
                ):
                    # Perform the gene wise quantification for a given library.
                    gene_quanti_per_lib_species_folder = (
                        self._pathcreator.gene_quanti_folders_by_species[sp][
                            "gene_quanti_per_lib_folder"
                        ]
                    )
                    gene_quanti_paths = [
                        self._pathcreator.gene_quanti_paths_by_species(
                            gene_quanti_per_lib_species_folder,
                            lib_name,
                            annotation_file,
                        )
                        for annotation_file in annotation_files
                    ]
                    # Check if all output files for this library exist - if so
                    # skip their creation
                    if not any(
                        [
                            self._file_needs_to_be_created(
                                gene_quanti_path, quiet=True
                            )
                            for gene_quanti_path in gene_quanti_paths
                        ]
                    ):
                        sys.stderr.write(
                            "The file(s) %s exist(s). Skipping their/its generation.\n"
                            % ", ".join(gene_quanti_paths)
                        )
                        continue
                    strand_specific = True
                    if self._args.non_strand_specific:
                        strand_specific = False

                    if not self._args.count_cross_aligned_reads:
                        crossmapped_reads = self._crossmapped_reads_by_lib[
                            lib_name
                        ]
                    else:
                        crossmapped_reads = None

                    gene_wise_quantification = GeneWiseQuantification(
                        references_by_species=references_by_species,
                        min_overlap=self._args.min_overlap,
                        read_region=self._args.read_region,
                        clip_length=self._args.clip_length,
                        norm_by_alignment_freq=norm_by_alignment_freq,
                        norm_by_overlap_freq=norm_by_overlap_freq,
                        allowed_features_str=self._args.allowed_features,
                        add_antisense=self._args.add_antisense,
                        antisense_only=self._args.antisense_only,
                        strand_specific=strand_specific,
                        unique_only=self._args.unique_only,
                        count_cross_aligned_reads=self._args.count_cross_aligned_reads,
                        crossmapped_reads=crossmapped_reads,
                    )

                    if norm_by_overlap_freq:
                        gene_wise_quantification.calc_overlaps_per_alignment(
                            read_alignment_path,
                            self._pathcreator.annotation_paths_by_species[sp],
                        )
                    for annotation_file, annotation_path in zip(
                        annotation_files,
                        self._pathcreator.annotation_paths_by_species[sp],
                    ):
                        output_path = (
                            self._pathcreator.gene_quanti_paths_by_species(
                                gene_quanti_per_lib_species_folder,
                                lib_name,
                                annotation_file,
                            )
                        )
                        jobs.append(
                            executor.submit(
                                gene_wise_quantification.quantify,
                                read_alignment_path,
                                annotation_path,
                                output_path,
                                self._args.pseudocounts,
                            )
                        )
        # Evaluate thread outcome
        self._check_job_completeness(jobs)

        for (
            sp,
            annotation_files,
        ) in self._pathcreator.annotation_files_by_species.items():
            self._gene_quanti_create_overview(
                sp,
                annotation_files,
                self._pathcreator.annotation_paths_by_species[sp],
                lib_names,
                norm_by_fragments
            )

    def _was_paired_end_alignment(self, lib_names):
        """Check if the mapping was done in paired- or single-end mode"""
        if len(lib_names) * 2 == len(self._pathcreator.get_read_files()):
            return True
        return False

    def _gene_quanti_create_overview(
        self, sp, annotation_files, annotation_paths, lib_names, norm_by_fragments
    ):
        """Create an overview table of all gene quantification for all libs."""
        strand_specific = True
        if self._args.non_strand_specific:
            strand_specific = False
        gene_wise_overview = GeneWiseOverview(
            allowed_features_str=self._args.allowed_features,
            add_antisense=self._args.add_antisense,
            antisense_only=self._args.antisense_only,
            strand_specific=strand_specific,
        )

        gene_quanti_per_lib_species_folder = (
            self._pathcreator.gene_quanti_folders_by_species[sp][
                "gene_quanti_per_lib_folder"
            ]
        )

        path_and_name_combos = {}
        for annotation_file, annotation_path in zip(
            annotation_files, annotation_paths
        ):
            path_and_name_combos[annotation_path] = []
            for lib in lib_names:
                path_and_name_combos[annotation_path].append(
                    [
                        lib,
                        self._pathcreator.gene_quanti_paths_by_species(
                            gene_quanti_per_lib_species_folder,
                            lib,
                            annotation_file,
                        ),
                    ]
                )

        if self._file_needs_to_be_created(
            self._pathcreator.gene_quanti_files_by_species[sp][
                "gene_wise_quanti_combined_path"
            ]
        ):
            gene_wise_overview.create_overview_raw_countings(
                path_and_name_combos,
                lib_names,
                self._pathcreator.gene_quanti_files_by_species[sp][
                    "gene_wise_quanti_combined_path"
                ],
            )
        if self._file_needs_to_be_created(
            self._pathcreator.gene_quanti_files_by_species[sp][
                "gene_wise_quanti_combined_rpkm_path"
            ]
        ):
            gene_wise_overview.create_overview_rpkm(
                path_and_name_combos,
                lib_names,
                self._pathcreator.gene_quanti_files_by_species[sp][
                    "gene_wise_quanti_combined_rpkm_path"
                ],
                self._libs_and_total_num_of_aligned_reads(
                    sp, self._args.normalize_cross_aligned_reads_included, norm_by_fragments
                ),
            )
        if self._file_needs_to_be_created(
            self._pathcreator.gene_quanti_files_by_species[sp][
                "gene_wise_quanti_combined_tnoar_path"
            ]
        ):
            gene_wise_overview.create_overview_norm_by_tnoar(
                path_and_name_combos,
                lib_names,
                self._pathcreator.gene_quanti_files_by_species[sp][
                    "gene_wise_quanti_combined_tnoar_path"
                ],
                self._libs_and_total_num_of_aligned_reads(
                    sp, self._args.normalize_cross_aligned_reads_included, norm_by_fragments
                ),
            )
        if self._file_needs_to_be_created(
            self._pathcreator.gene_quanti_files_by_species[sp][
                "gene_wise_quanti_combined_tpm_path"
            ]
        ):
            gene_wise_overview.create_overview_tpm(
                self._pathcreator.gene_quanti_files_by_species[sp][
                    "gene_wise_quanti_combined_path"
                ],
                self._pathcreator.gene_quanti_files_by_species[sp][
                    "gene_wise_quanti_combined_tpm_path"
                ],
            )

    def _libs_and_total_num_of_aligned_reads(
        self, sp, normalize_cross_aligned_reads_included=False, norm_by_fragments=False,
    ):
        """Read the total number of reads per library for a selected species."""
        if norm_by_fragments:
            reads_or_fragments = "fragments"
            alignment_stats_path = self._pathcreator.fragment_alignments_stats_path
        else:
            reads_or_fragments = "reads"
            alignment_stats_path = self._pathcreator.read_alignments_stats_path
        with open(
            alignment_stats_path
        ) as read_aligner_stats_fh:
            read_aligner_stats = json.loads(read_aligner_stats_fh.read())
        if normalize_cross_aligned_reads_included:
            return dict(
                [
                    (lib, values["species_stats"][sp][f"no_of_aligned_{reads_or_fragments}"])
                    for lib, values in read_aligner_stats.items()
                ]
            )

        else:
            return dict(
                [
                    (
                        lib,
                        values["species_stats"][sp][f"no_of_aligned_{reads_or_fragments}"]
                        - values["species_stats"][sp][
                            f"no_of_cross_aligned_{reads_or_fragments}"
                        ],
                    )
                    for lib, values in read_aligner_stats.items()
                ]
            )

    def _libs_and_total_num_of_uniquely_aligned_reads(self):
        """Read the total number of reads per library."""
        with open(
            self._pathcreator.read_alignments_stats_path
        ) as read_aligner_stats_fh:
            read_aligner_stats = json.loads(read_aligner_stats_fh.read())
        return dict(
            [
                (lib, values["stats_total"]["no_of_uniquely_aligned_reads"])
                for lib, values in read_aligner_stats.items()
            ]
        )

    def compare_with_deseq(self):
        """Manage the pairwise expression comparison with DESeq."""
        project_creator = ProjectCreator()
        project_creator.create_subfolders(
            self._pathcreator.required_deseq_folders()
        )
        arg_libs = [
            self._pathcreator._clean_file_name(lib)
            for lib in self._args.libs.split(",")
        ]
        conditions = self._args.conditions.split(",")
        replicates = self._args.replicates.split(",")
        libs_by_species = self._get_libs_by_species(self._args.libs_by_species)
        size_factor = self._args.size_factor
        self._check_deseq_args(arg_libs, conditions)
        for sp in self._species_folder_prefixes_and_display_names.keys():
            deseq_runner = DESeqRunner(
                sp,
                arg_libs,
                conditions,
                replicates,
                libs_by_species,
                size_factor,
                self._pathcreator.deseq_folders_by_species[sp][
                    "deseq_raw_folder"
                ],
                self._pathcreator.deseq_folders_by_species[sp][
                    "deseq_extended_folder"
                ],
                self._pathcreator.deseq_files_by_species[sp][
                    "deseq_script_path"
                ],
                self._pathcreator.deseq_files_by_species[sp][
                    "deseq_pca_heatmap_path"
                ],
                self._pathcreator.gene_quanti_files_by_species[sp][
                    "gene_wise_quanti_combined_path"
                ],
                self._pathcreator.deseq_files_by_species[sp][
                    "deseq_tmp_session_info_script"
                ],
                self._pathcreator.deseq_files_by_species[sp][
                    "deseq_session_info"
                ],
                self._args.fc_shrinkage_off,
                self._args.cooks_cutoff_off,
            )
            deseq_runner.create_deseq_script_file()
            deseq_runner.write_session_info_file()
            deseq_runner.run_deseq()
            deseq_runner.merge_counting_files_with_results()

    def _get_libs_by_species(self, arg_libs_by_species):
        libs_by_species = {}
        for sp_and_libs in arg_libs_by_species:
            sp, libs_string = sp_and_libs.split("=")
            libs = libs_string.split(",")
            clean_libs = [
                self._pathcreator._clean_file_name(lib) for lib in libs
            ]
            libs_by_species[sp] = clean_libs
        return libs_by_species

    def _check_deseq_args(self, arg_libs, conditions):
        """Test if the given arguments are sufficient."""
        if len(arg_libs) != len(conditions):
            self._write_err_msg_and_quit(
                "Error - The read library file list and condition list must "
                "have the same number of elements. You entered \n%s "
                "(= %s elements)\nand \n%s (= %s elements).\n"
                % (
                    self._args.libs,
                    len(arg_libs),
                    self._args.conditions,
                    len(conditions),
                )
            )
        raw_stat_data_reader = RawStatDataReader()
        alignment_stats = [
            raw_stat_data_reader.read(
                self._pathcreator.read_alignments_stats_path
            )
        ]
        lib_names = list(alignment_stats[0].keys())
        if len(lib_names) != len(arg_libs):
            self._write_err_msg_and_quit(
                "The number of read libraries is lower or higher than "
                "expected. The following read libs are available: %s\nThe "
                'following read list string is suggested: "%s"\n'
                % (", ".join(lib_names), ",".join(lib_names))
            )
        for lib in lib_names:
            if lib not in arg_libs:
                self._write_err_msg_and_quit(
                    'The library "%s" is not present in your list of '
                    "libraries. Please add it.\n" % (lib)
                )

    def _write_err_msg_and_quit(self, msg):
        """Write error message and close the program gracefully."""
        sys.stderr.write(msg)
        sys.exit(1)

    def viz_align(self):
        """Generate plots based on the read processing and mapping"""
        # Create the output folders
        viz_align_folders = self._pathcreator.required_viz_align_folders()
        project_creator = ProjectCreator()
        project_creator.create_subfolders(viz_align_folders)
        align_viz = AlignViz(
            self._pathcreator.get_lib_names_single_end()
            if not self._args.paired_end
            else self._pathcreator.get_lib_names_paired_end(),
            self._pathcreator.read_processing_stats_path,
            self._pathcreator.read_alignments_stats_path,
            self._pathcreator.read_alignment_stats_table_path,
            self._pathcreator.viz_align_aligned_reads_by_species_paths,
            self._species_folder_prefixes_and_display_names,
        )
        align_viz.read_stat_files()
        align_viz.plot_input_read_length(
            self._pathcreator.viz_align_input_read_length_plot_path
        )
        align_viz.plot_processed_read_length(
            self._pathcreator.viz_align_processed_reads_length_plot_path
        )
        align_viz.plot_total_number_of_aligned_reads(
            self._pathcreator.viz_align_all_folder
        )
        align_viz.plot_species_exclusive_reads_for_each_species()

    def viz_gene_quanti(self):
        """Generate plots based on the gene-wise read countings"""
        # Create output folders for each species
        project_creator = ProjectCreator()
        project_creator.create_subfolders(
            self._pathcreator.required_viz_gene_quanti_folders()
        )

        for sp in self._species_folder_prefixes:
            # Set output folder and files paths for each species
            gene_wise_quanti_combined_path = (
                self._pathcreator.gene_quanti_files_by_species[sp][
                    "gene_wise_quanti_combined_path"
                ]
            )
            viz_gene_quanti_scatter_plot_path = (
                self._pathcreator.viz_gene_quanti_files_by_species[sp][
                    "viz_gene_quanti_scatter_plot_path"
                ]
            )
            rna_classes_plot_path = (
                self._pathcreator.viz_gene_quanti_files_by_species[sp][
                    "rna_classes_plot_path"
                ]
            )

            # Create plots
            gene_quanti_viz = GeneQuantiViz(
                gene_wise_quanti_combined_path,
                self._pathcreator.get_lib_names_single_end()
                if not self._args.paired_end
                else self._pathcreator.get_lib_names_paired_end(),
            )
            gene_quanti_viz.parse_input_table()
            gene_quanti_viz.plot_correlations(viz_gene_quanti_scatter_plot_path)
            gene_quanti_viz.plot_annotation_class_quantification(
                rna_classes_plot_path
            )

    def viz_deseq(self):
        """Generate plots based on the DESeq analysis"""
        # Create output folders for each species
        project_creator = ProjectCreator()
        project_creator.create_subfolders(
            self._pathcreator.required_viz_deseq_folders()
        )

        for sp in self._species_folder_prefixes:
            # Set output folder and files paths for each species
            deseq_path_template = (
                self._pathcreator.deseq_folders_by_species[sp][
                    "deseq_raw_folder"
                ]
                + "/deseq_comp_"
            )

            deseq_viz = DESeqViz(
                self._pathcreator.deseq_files_by_species[sp][
                    "deseq_script_path"
                ],
                deseq_path_template,
                max_pvalue=self._args.max_pvalue,
            )
            deseq_viz.create_scatter_plots(
                self._pathcreator.viz_deseq_files_by_species[sp][
                    "viz_deseq_scatter_plot_path"
                ]
            )
            deseq_viz.create_volcano_plots(
                self._pathcreator.viz_deseq_files_by_species[sp][
                    "viz_deseq_volcano_plot_path"
                ],
                self._pathcreator.viz_deseq_files_by_species[sp][
                    "viz_deseq_volcano_plot_adj_path"
                ],
            )
