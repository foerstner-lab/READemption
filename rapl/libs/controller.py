import concurrent.futures
import os
import sys
sys.path.append(".")
from libs.fasta import FastaParser
from libs.coveragecreator import CoverageCreator
from libs.parameterlog import ParameterLogger
from libs.parameters import Parameters
from libs.paths import Paths
from libs.projectcreator import ProjectCreator
from libs.readclipper import ReadClipper
from libs.readmapper import ReadMapper
from libs.readmapperstats import ReadMapperStats, ReadMapperStatsReader
from libs.seqsizefilter import SeqSizeFilter
from libs.annotationoverlap import AnnotationOverlap
from libs.annotationoverview import AnnotationOverview
from libs.sambamconverter import SamToBamConverter

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
        self.paths.set_read_files_dep_file_lists(read_file_names)
        self.paths.set_genome_paths(genome_file_names)
        # TODO
        # - generate input stats before mapping
        # - read tracing after mapping
        self._prepare_reads()
        self._map_reads()
        self._sam_to_bam()
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
            self.paths.read_mapping_result_sam_paths,
            self.paths.unmapped_reads_paths,
            int(self.args.threads),
            int(self.args.segemehl_accuracy),
            int(self.args.segemehl_evalue))

    def _sam_to_bam(self):
        sam_to_bam_converter = SamToBamConverter(
            samtools_bin=self.args.samtools_bin)
        for sam_path, bam_prefix_path in zip(
            self.paths.read_mapping_result_sam_paths,
            self.paths.read_mapping_result_bam_prefixes_paths):
            sam_to_bam_converter.sam_to_bam(sam_path, bam_prefix_path)

    def _generate_read_mapping_stats(self, read_file_names):
        ref_ids_to_file_name = self._ref_ids_to_file_name(
            self.paths.genome_file_paths)
        read_mapper_stats = ReadMapperStats(
            samtools_bin=self.args.samtools_bin)
        read_mapper_stats.count_raw_reads(
            read_file_names, self.paths.read_file_paths)
        read_mapper_stats.count_long_enough_clipped_reads(
            read_file_names, self.paths.clipped_read_file_long_enough_paths)
        read_mapper_stats.count_too_small_clipped_reads(
            read_file_names, self.paths.clipped_read_file_too_short_paths)
        read_mapper_stats.count_mappings(
            read_file_names, self.paths.read_mapping_result_bam_paths)
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

    def create_coverage_files(self):
        """Create coverage files based on the Segemehl mappings."""
        read_file_names = self.paths._get_read_file_names()
        self.paths.set_read_files_dep_file_lists(read_file_names)
        read_mapper_stat_reader = ReadMapperStatsReader()
        read_mapping_stats = read_mapper_stat_reader.read_stat_file(
            self.paths.read_mapping_stat_file)
        min_read_mapping_counting = read_mapper_stat_reader.min_read_countings(
            self.paths.read_mapping_stat_file)
        for read_file_name, bam_file_path in zip(
            read_file_names, self.paths.read_mapping_result_bam_paths):
            self._create_coverage_files_for_lib(
                read_file_name, bam_file_path, read_mapping_stats,
                min_read_mapping_counting)

        # TODO
        # # Run the generation of coverage in parallel
        # threads = []
        # with concurrent.futures.ThreadPoolExecutor(
        #     max_workers=self.parameters.python_number_of_threads) as executor:
        #     for read_file_name, bam_file_path in zip(
        #         read_file_names, self.paths.read_mapping_result_bam_paths):
        #         threads.append(executor.submit(
        #                 self._create_coverage_files_for_lib,
        #                 read_file_name, bam_file_path, read_mapping_stats,
        #                 min_read_mapping_counting))
        # Evaluate thread outcome
        # self._check_thread_completeness(threads)

    def _create_coverage_files_for_lib(
        self, read_file_name, bam_file_path, read_mapping_stats,
        min_read_mapping_counting):
        coverage_creator = CoverageCreator(
            samtools_bin=self.args.samtools_bin)
        read_count_splitting = True
        if self.args.skip_read_count_splitting:
            read_count_splitting = False
        coverage_creator.init_coverage_lists(bam_file_path)
        coverage_creator.count_coverage(
            bam_file_path, read_count_splitting=read_count_splitting,
            uniqueley_mapped_only=self.args.unique_only)
        # Raw countings
        coverage_creator.write_to_files(
            "%s/%s" % (self.paths.coverage_folder, read_file_name),
            read_file_name)
        total_number_of_mapped_reads = read_mapping_stats[read_file_name][
            "total_number_of_mapped_reads"]
        # Read normalized countings - multiplied by min read counting
        factor = (min_read_mapping_counting / total_number_of_mapped_reads)
        coverage_creator.write_to_files(
            "%s/%s-div_by_%.1f_multi_by_%.1f" % (
                self.paths.coverage_folder_norm_reads, read_file_name,
                total_number_of_mapped_reads, min_read_mapping_counting),
            read_file_name, factor=factor)
        # Read normalized countings - multiplied by 1M
        factor = (1000000 / total_number_of_mapped_reads)
        coverage_creator.write_to_files(
            "%s/%s-div_by_%.1f_multi_by_1M" % (
                self.paths.coverage_folder_norm_reads_mil, read_file_name,
                total_number_of_mapped_reads), read_file_name, factor=factor)

    def search_annotation_overlaps(self):
        """Search for overlaps of reads and annotations."""
        annotation_overlap = AnnotationOverlap()
        read_file_names = self.paths._get_read_file_names()
        annotation_files = self.paths._get_annotation_file_names()
        self.paths.set_annotation_paths(annotation_files)
        self.paths.set_read_files_dep_file_lists(read_file_names)
        threads = []

        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.args.threads) as executor:
            for read_file_name, read_mapping_path in zip(
                read_file_names, self.paths.read_mapping_result_bam_paths):
                for annotation_file, annotation_file_path in zip(
                    annotation_files, self.paths.annotation_file_paths):
                    annotation_hit_file_path = (
                        self.paths.annotation_hit_file_path(
                            read_file_name, annotation_file))
                    threads.append(executor.submit(
                            annotation_overlap.find_overlaps, read_mapping_path,
                            annotation_file_path, annotation_hit_file_path))
        self._check_thread_completeness(threads)

    def _check_thread_completeness(self, threads):
        """Check the completness of each thread in a list"""
        for thread in concurrent.futures.as_completed(threads):
            if thread.exception():
                raise(thread.exception())

    def create_annotation_overview(self):
        self.annotation_overview = AnnotationOverview(
            min_overlap=self.args.min_overlap)
        read_file_names = self.paths._get_read_file_names()
        annotation_files = self.paths._get_annotation_file_names()
        self.paths.set_annotation_paths(annotation_files)
        self.paths.set_read_files_dep_file_lists(read_file_names)
        self._count_annotation_hits(read_file_names, annotation_files)
        self._write_anno_overview_files(read_file_names, annotation_files)

    def _count_annotation_hits(self, read_file_names, annotation_files):
        for read_file_name, read_mapping_path in zip(
            read_file_names, self.paths.read_mapping_result_bam_paths):
            # Get the number of read mappings for each read - this
            # needs only be done once for each read file
            self.annotation_overview.get_read_mapping_freq(read_mapping_path)
            self.annotation_overview.add_mapping_counting_table(read_file_name)
            for annotation_file, annotation_file_path in zip(
                annotation_files, self.paths.annotation_file_paths):
                # Get the number of overlaps each read mapping
                # has. Sense and antisense are considered.
                annotation_hit_file_path = self.paths.annotation_hit_file_path(
                    read_file_name, annotation_file)
                self.annotation_overview.get_mapping_overlap_freq(
                    read_file_name, annotation_hit_file_path)
            for annotation_file, annotation_file_path in zip(
                annotation_files, self.paths.annotation_file_paths):
                self.annotation_overview.init_counting_table(
                    read_file_name, annotation_file, annotation_file_path)
                annotation_hit_file_path = self.paths.annotation_hit_file_path(
                    read_file_name, annotation_file)
                self.annotation_overview.count_overlapping_reads_per_gene(
                    read_file_name, annotation_file, annotation_hit_file_path)
            # Once all annotation hit files are processed clean the
            # associated dictionaries
            self.annotation_overview.clean_dicts()

    def _write_anno_overview_files(self, read_file_names, annotation_files):
        for annotation_file, annotation_file_path in zip(
                annotation_files, self.paths.annotation_file_paths):
            annotation_hit_overview_sense_file_path = (
                self.paths.annotation_hit_overview_file_path(
                    annotation_file, "sense"))
            annotation_hit_overview_antisense_file_path = (
                self.paths.annotation_hit_overview_file_path(
                    annotation_file, "antisense"))
            # Raw countings
            self.annotation_overview.write_overview_tables(
                annotation_file, annotation_file_path, read_file_names,
                annotation_hit_overview_sense_file_path,
                annotation_hit_overview_antisense_file_path)
            # Mapped reads normalized countings
            # read_mapper_stat_reader = ReadMapperStatsReader()
            # read_mapping_stats = read_mapper_stat_reader.read_stat_file(
            #     self.paths.read_mapping_stat_file)
            # total_numbers_of_mapped_reads = [
            #     read_mapping_stats[read_file_name][
            #         "total_number_of_mapped_reads"]
            #     for read_file_names in read_file_names]
            # print(total_number_of_mapped_reads)
            # self.annotation_overview.write_overview_tables(
            #     annotation_file, annotation_file_path, read_file_names,
            #     annotation_hit_overview_sense_file_path,
            #     annotation_hit_overview_antisense_file_path,
            #     total_numbers_of_mapped_reads=total_numbers_of_mapped_reads)




