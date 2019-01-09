from reademptionlib.deseq2 import DESeq2Runner
from reademptionlib.helpers import Helpers
from reademptionlib.paths import Paths
from reademptionlib.rawstatdata import RawStatDataReader


class DESeq2Controller(object):
    
    def __init__(self, args):
        self._args = args
        self._paths = Paths(args)
        self._helpers = Helpers(args)

    def compare_with_deseq2(self):
        """Manage the pairwise expression comparison with DESeq2."""
        self._helpers.test_folder_existance(
            self._paths.required_deseq2_folders())
        arg_libs = [self._paths._clean_file_name(lib) for lib in
                    self._args.libs.split(",")]
        conditions = self._args.conditions.split(",")
        self._check_deseq2_args(arg_libs, conditions)
        deseq2_runner = DESeq2Runner(
            arg_libs, conditions, self._paths.deseq2_raw_folder,
            self._paths.deseq2_extended_folder, self._paths.deseq2_script_path,
            self._paths.deseq2_pca_heatmap_path,
            self._paths.gene_wise_quanti_combined_path,
            self._paths.deseq2_tmp_session_info_script,
            self._paths.deseq2_session_info,
            self._args.cooks_cutoff_off)
        deseq2_runner.create_deseq2_script_file()
        deseq2_runner.write_session_info_file()
        deseq2_runner.run_deseq2()
        deseq2_runner.merge_counting_files_with_results()
        deseq2_runner.create_final_output_files()

    def _check_deseq2_args(self, arg_libs, conditions):
        """Test if the given arguments are sufficient."""
        if len(arg_libs) != len(conditions):
            self._helpers.write_err_msg_and_quit(
                "Error - The read library file list and condition list must "
                "have the same number of elements. You entered \n%s "
                "(= %s elements)\nand \n%s (= %s elements).\n" % (
                    self._args.libs, len(arg_libs), self._args.conditions,
                    len(conditions)))
        raw_stat_data_reader = RawStatDataReader()
        alignment_stats = [raw_stat_data_reader.read(
            self._paths.read_alignments_stats_path)]
        lib_names = list(alignment_stats[0].keys())
        if len(lib_names) != len(arg_libs):
            self._helpers.write_err_msg_and_quit(
                "The number of read libraries is lower or higher than "
                "expected. The following read libs are available: %s\nThe "
                "following read list string is suggested: \"%s\"\n" % (
                    ", ".join(lib_names), ",".join(lib_names)))
        for lib in lib_names:
            if lib not in arg_libs:
                self._helpers.write_err_msg_and_quit(
                    "The library \"%s\" is not present in your list of "
                    "libraries. Please add it.\n" % (lib))
