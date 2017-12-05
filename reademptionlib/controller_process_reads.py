import concurrent.futures
from reademptionlib.cutadapt import Cutadapt
from reademptionlib.helpers import Helpers
from reademptionlib.paths import Paths
# from reademptionlib.readprocessor import ReadProcessor


class ReadProcessingController():

    def __init__(self, args):
        """Create an instance."""
        self._args = args
        self._paths = Paths(args)
        self._helpers = Helpers(args)

    def process_reads(self):
        if not self._args.paired_end:
            self._read_files = self._paths.get_read_files()
            self._lib_names = self._paths.get_lib_names_single_end()
            self._paths.set_read_files_dep_file_lists_single_end(
                self._read_files, self._lib_names)
            self._prepare_reads_se_cutadapt()
        else:
            self._prepare_reads_pe_cutadapt()

    def _prepare_reads_se_cutadapt(self):
        read_files_and_jobs = {}
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=self._args.processes) as executor:
            for lib_name, read_path, processed_read_path in zip(
                    self._lib_names, self._paths.read_paths,
                    self._paths.processed_read_paths):
                if not self._helpers.file_needs_to_be_created(
                        processed_read_path):
                    continue
                cutadapt = Cutadapt(
                    self._args, self._args.cutadapt_options,
                    self._args.cutadapt_bin)
                read_files_and_jobs[lib_name] = executor.submit(
                    cutadapt.run_cutadapt_se, read_path,
                    self._paths.read_processing_base_folder, lib_name)
        # self._paths.gzip_processed_reads()
        # self._helpers.check_job_completeness(read_files_and_jobs.values())

    def _prepare_reads_pe_cutadapt(self):
        read_files_and_jobs = {}
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=self._args.processes) as executor:
            for lib_name, read_path_pair, processed_read_path_pair in zip(
                    self._lib_names, self._paths.read_path_pairs,
                    self._paths.processed_read_path_pairs):
                if not self._helpers.file_needs_to_be_created(
                        processed_read_path_pair):
                    continue
                cutadapt = Cutadapt(
                    self._args.cutadapt_options, self._args.cutadapt_bin)
                read_files_and_jobs[lib_name] = executor.submit(
                    cutadapt.run_cutadapt_pe, read_path_pair,
                    self._paths.processed_reads_folder, lib_name)
        #self._helpers.check_job_completeness(read_files_and_jobs.values())
        

    # def _prepare_reads_single_end(self):
    #     """Manage the prepartion of reads before the actual mappings."""
    #     read_files_and_jobs = {}
    #     with concurrent.futures.ProcessPoolExecutor(
    #             max_workers=self._args.processes) as executor:
    #         for lib_name, read_path, processed_read_path in zip(
    #                 self._lib_names, self._paths.read_paths,
    #                 self._paths.processed_read_paths):
    #             if not self._helpers.file_needs_to_be_created(
    #                     processed_read_path):
    #                 continue
    #             read_processor = ReadProcessor(
    #                 poly_a_clipping=self._args.poly_a_clipping,
    #                 min_read_length=self._args.min_read_length,
    #                 min_phred_score=self._args.min_phred_score,
    #                 adapter=self._args.adapter,
    #                 reverse_complement=self._args.reverse_complement)
    #             read_files_and_jobs[lib_name] = executor.submit(
    #                 read_processor.process_single_end, read_path,
    #                 processed_read_path)
    #     self._evaluet_job_and_generate_stat_file(read_files_and_jobs)


    # def _prepare_reads_paired_end(self):
    #     read_files_and_jobs = {}
    #     with concurrent.futures.ProcessPoolExecutor(
    #             max_workers=self._args.processes) as executor:
    #         for lib_name, read_path_pair, processed_read_path_pair in zip(
    #             self._lib_names, self._paths.read_path_pairs,
    #                 self._paths.processed_read_path_pairs):
    #             for processed_read_path in processed_read_path_pair:
    #                 if not self._helpers.file_needs_to_be_created(
    #                         processed_read_path):
    #                     continue
    #                 read_processor = ReadProcessor(
    #                     poly_a_clipping=False,
    #                     min_read_length=self._args.min_read_length,
    #                     min_phred_score=self._args.min_phred_score,
    #                     adapter=self._args.adapter)
    #                 read_files_and_jobs[lib_name] = executor.submit(
    #                     read_processor.process_paired_end, read_path_pair,
    #                     processed_read_path_pair)
    #     self._evaluet_job_and_generate_stat_file(read_files_and_jobs)

