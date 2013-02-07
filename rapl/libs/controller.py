import concurrent.futures
import os
import sys
sys.path.append(".")
from libs.fasta import FastaParser
from libs.coveragecreator import CoverageCreator
from libs.parameterlog import ParameterLogger
from libs.paths import Paths
from libs.projectcreator import ProjectCreator
from libs.readprocessor import ReadProcessor
from libs.readmapper import ReadMapper
from libs.readmapperstats import ReadMapperStats, ReadMapperStatsReader
from libs.seqsizefilter import SeqSizeFilter
from libs.sambamconverter import SamToBamConverter
from libs.genewisequanti import GeneWiseQuantification, GeneWiseOverview
from libs.rawstatdata import RawStatDataWriter, RawStatDataReader

class Controller(object):

    def __init__(self, args):
        """Create an instance."""
        self.args = args
        self.paths = Paths(args.project_path)

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
        self.read_file_names = self.paths._get_read_file_names()
        genome_file_names = self.paths._get_genome_file_names()
        self.paths.set_read_files_dep_file_lists(self.read_file_names)
        self.paths.set_genome_paths(genome_file_names)
        self._prepare_reads()
        self._map_reads()
        self._sam_to_bam()
        self._generate_read_mapping_stats()
        self._write_mapping_stat_table()

    def _prepare_reads(self):
        raw_stat_data_writer = RawStatDataWriter(pretty=True)
        combined_stats = {}
        read_processor = ReadProcessor(min_read_length=self.args.min_read_length)
        for read_file, read_file_path, processed_read_file_path in zip(
                self.read_file_names,
                self.paths.read_file_paths,
                self.paths.processed_read_file_paths):
            processing_stats = read_processor.process(
                read_file_path, processed_read_file_path)
            combined_stats[read_file] = processing_stats
        raw_stat_data_writer.write(
            combined_stats, self.paths.read_processing_stats_path)

    def _map_reads(self):
        read_mapper = ReadMapper(segemehl_bin=self.args.segemehl_bin)
        read_mapper.build_index(
            self.paths.genome_file_paths, self.paths.index_file_path)
        read_mapper.run_mappings(
            self.paths.processed_read_file_paths,
            self.paths.genome_file_paths, self.paths.index_file_path,
            self.paths.read_mapping_result_sam_paths,
            self.paths.unmapped_reads_paths,
            int(self.args.threads),
            int(self.args.segemehl_accuracy),
            int(self.args.segemehl_evalue))

    def _sam_to_bam(self):
        sam_to_bam_converter = SamToBamConverter()
        for sam_path, bam_prefix_path in zip(
            self.paths.read_mapping_result_sam_paths,
            self.paths.read_mapping_result_bam_prefixes_paths):
            sam_to_bam_converter.sam_to_bam(sam_path, bam_prefix_path)

    def _generate_read_mapping_stats(self):
        raw_stat_data_writer = RawStatDataWriter(pretty=True)
        combined_stats = {}
        for (read_file_name, read_mapping_result_bam_path,
             unmapped_reads_path) in zip(
                self.read_file_names,
                self.paths.read_mapping_result_bam_paths,
                self.paths.unmapped_reads_paths):
            read_mapping_stats = ReadMapperStats()
            combined_stats[read_file_name] = read_mapping_stats.count(
                read_mapping_result_bam_path,
                unmapped_reads_path)
        raw_stat_data_writer.write(
            combined_stats, self.paths.read_mapping_stats_path)

    def _write_mapping_stat_table(self):
        raw_stat_data_reader = RawStatDataReader()
        read_processing_stats = raw_stat_data_reader.read(
            self.paths.read_processing_stats_path)
        mapping_stats = raw_stat_data_reader.read(
            self.paths.read_mapping_stats_path)
        table = []
        table.append(["Lib"] + self.read_file_names)
        ref_ids = sorted(list(list(mapping_stats.values())[0][
            "countings_per_reference"].keys()))
        table = [
            ["Library read file"] + self.read_file_names,
            ["No. of input reads"] +
            self._get_read_process_numbers(
                read_processing_stats, "total_no_of_reads"),
            [ "No. of reads - PolyA detected and removed"] +
            self._get_read_process_numbers(
                read_processing_stats, "polya_removed"),
            ["No. of reads - Single 3' A removed"] +
            self._get_read_process_numbers(
                read_processing_stats, "single_a_removed"),
            ["No. of reads - Unmodified"] + self._get_read_process_numbers(
            read_processing_stats, "unmodified"),
            ["No. of reads - Removed as too short"] +
            self._get_read_process_numbers(read_processing_stats, "too_short"),
            ["No. of reads - Long enough and used for mapping"] +
            self._get_read_process_numbers(
                read_processing_stats, "long_enough"),
            ["Total no. of mapped read"] + [
                round(num) for num in self._total_mapping_stat_numbers(
                mapping_stats, "no_of_mapped_reads")],
            ["Total no. of unmapped reads"] + [str(mapping_stats[read_file][
                "no_of_unmapped_reads"]) for read_file in self.read_file_names],
            ["Total no. of uniquely mapped reads"] +
            self._total_mapping_stat_numbers(
                mapping_stats, "no_of_uniquely_mapped_reads"),
            ["Total no. of mappings"] + self._total_mapping_stat_numbers(
                mapping_stats, "no_of_mappings"),
            ["Percentag of mapped reads (compared to total input reads)"]  + [
                round(float(mapped_reads)/float(total_reads) * 100, 2)
                for mapped_reads, total_reads in
                zip(self._total_mapping_stat_numbers(
                    mapping_stats, "no_of_mapped_reads"),
                    self._get_read_process_numbers(
                        read_processing_stats, "total_no_of_reads"))],
            ["Percentag of uniquely mapped reads (in relation to all mapped reads)"]  + [
                round(float(uniquely_mapped_reads)/float(mapped_reads) * 100, 2)
                for uniquely_mapped_reads, mapped_reads  in zip(
                        self._total_mapping_stat_numbers(mapping_stats, "no_of_uniquely_mapped_reads"),
                        self._total_mapping_stat_numbers(mapping_stats, "no_of_mapped_reads"))]
        ]
        for ref_id in ref_ids:
            table.append(
                ["%s - No. of mapped reads" % ref_id] +
                self._mapping_number_per_ref_seq(
                    mapping_stats, ref_id, "no_of_mapped_reads"))
            table.append(
                ["%s - No. of uniquely mapped reads" % ref_id] +
                self._mapping_number_per_ref_seq(
                    mapping_stats, ref_id, "no_of_uniquely_mapped_reads"))
            table.append(
                ["%s - No. of mapping" % ref_id] +
                self._mapping_number_per_ref_seq(
                    mapping_stats, ref_id, "no_of_mappings"))
        table_fh = open(self.paths.read_mapping_stats_table_path, "w")
        table_fh.write("\n".join(["\t".join([str(cell) for cell in row]) for row in table]))
        table_fh.close()

    def _mapping_number_per_ref_seq(self, mapping_stats, ref_id, attribute):
        return([mapping_stats[read_file]["countings_per_reference"][
            ref_id]["no_of_mapped_reads"] for read_file in self.read_file_names])

    def _total_mapping_stat_numbers(self, mapping_stats, attribute):
        return([mapping_stats[read_file]["countings_total"][attribute]
                for read_file in self.read_file_names])

    def _get_read_process_numbers(
            self, read_processing_stats, attribute):
        return([read_processing_stats[read_file][attribute]
                for read_file in self.read_file_names])

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
        # Run the generation of coverage in parallel
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self.args.threads) as executor:
            for read_file_name, bam_file_path in zip(
                read_file_names, self.paths.read_mapping_result_bam_paths):
                jobs.append(executor.submit(
                        self._create_coverage_files_for_lib,
                        read_file_name, bam_file_path, read_mapping_stats,
                        min_read_mapping_counting))
        # Evaluate thread outcome
        self._check_job_completeness(jobs)

    def _create_coverage_files_for_lib(
        self, read_file_name, bam_file_path, read_mapping_stats,
        min_read_mapping_counting):
        coverage_creator = CoverageCreator()
        read_count_splitting = True
        if self.args.skip_read_count_splitting is True:
            read_count_splitting = False
        coverage_creator.init_coverage_lists(bam_file_path)
        coverage_creator.count_coverage(
            bam_file_path, read_count_splitting=read_count_splitting,
            uniqueley_mapped_only=self.args.unique_only,
            first_base_only=self.args.first_base_only)
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

    def _check_job_completeness(self, jobs):
        """Check the completness of each job in a list"""
        for job in concurrent.futures.as_completed(jobs):
            if job.exception():
                raise(job.exception())

    def quantify_gene_wise(self):
        norm_by_mapping_freq = True
        norm_by_overlap_freq = True
        if self.args.skip_norm_by_mapping_freq:
            norm_by_mapping_freq = False
        if self.args.skip_norm_by_overlap_freq:
            norm_by_overlap_freq = False
        read_file_names = self.paths._get_read_file_names()
        annotation_files = self.paths._get_annotation_file_names()
        self.paths.set_annotation_paths(annotation_files)
        self.paths.set_read_files_dep_file_lists(read_file_names)
        for read_file_name, read_mapping_path in zip(
                read_file_names, self.paths.read_mapping_result_bam_paths):
            gene_wise_quantification = GeneWiseQuantification(
                min_overlap=self.args.min_overlap,
                norm_by_mapping_freq=norm_by_mapping_freq,
                norm_by_overlap_freq=norm_by_overlap_freq)
            gene_wise_quantification.calc_overlaps_per_mapping(
                read_mapping_path, self.paths.annotation_file_paths)
            for  annotation_file, annotation_file_path in zip(
                    annotation_files, self.paths.annotation_file_paths):
                gene_wise_quantification.quantify(
                    read_mapping_path, annotation_file_path,
                    self.paths.gene_quanti_path(
                        read_file_name, annotation_file))
        self._gene_quanti_create_overview(
            annotation_files, self.paths.annotation_file_paths, read_file_names)

    def _gene_quanti_create_overview(
            self, annotation_files, annotation_file_paths, read_file_names):
        gene_wise_overview = GeneWiseOverview()
        path_and_name_combos = {}
        for annotation_file, annotation_file_path in zip(
                annotation_files, annotation_file_paths):
            path_and_name_combos[annotation_file_path] = []
            for read_file_name in read_file_names:
                path_and_name_combos[annotation_file_path].append(
                    [read_file_name, self.paths.gene_quanti_path(
                        read_file_name, annotation_file)])
        gene_wise_overview.create_overview(
            path_and_name_combos, read_file_names,
            self.paths.gene_wise_quanti_combined_path)
