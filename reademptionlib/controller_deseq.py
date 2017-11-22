from reademptionlib.deseq import DESeqRunner
from reademptionlib.helpers import Helpers
from reademptionlib.paths import Paths
from reademptionlib.rawstatdata import RawStatDataReader


class RunDeseq(object):
    def __init__(self, args):
        self._args = args
        self._paths = Paths(args)
        self._helpers = Helpers(args)

    def compare_with_deseq(self):
        """Manage the pairwise expression comparison with DESeq."""
        self._helpers.test_folder_existance(
            self._paths.required_deseq_folders())
        arg_libs = [self._paths._clean_file_name(lib) for lib in
                    self._args.libs.split(",")]
        conditions = self._args.conditions.split(",")
        self._check_deseq_args(arg_libs, conditions)
        deseq_runner = DESeqRunner(
            arg_libs, conditions, self._paths.deseq_raw_folder,
            self._paths.deseq_extended_folder, self._paths.deseq_script_path,
            self._paths.deseq_pca_heatmap_path,
            self._paths.gene_wise_quanti_combined_path,
            self._paths.deseq_tmp_session_info_script,
            self._paths.deseq_session_info,
            self._args.cooks_cutoff_off)
        deseq_runner.create_deseq_script_file()
        deseq_runner.write_session_info_file()
        deseq_runner.run_deseq()
        deseq_runner.merge_counting_files_with_results()
        self._viz_deseq()
        deseq_runner.create_final_output_files()

    def _check_deseq_args(self, arg_libs, conditions):
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

    def _viz_deseq(self):
        """Generate plots based on the DESeq analysis"""
        from reademptionlib.vizdeseq import DESeqViz
        conditions = list(set(self._args.conditions.split(',')))
        comparison_path_template_1 = "{}/deseq_comp_{}_vs_{}_with_annotation_and_countings.csv".format(
            self._paths.deseq_extended_folder, conditions[0], conditions[1])
        comparison_1 = "{}_vs_{}".format(conditions[0], conditions[1])
        deseq_viz = DESeqViz(
            comparison_path_template_1, self._paths.viz_deseq_base_folder,
            self._args.padj_cutoff, comparison_1, self._args.alpha,
            self._args.color_sig, self._args.color_non_sig, self._args.shape,
            self._args.glyph_size)
        deseq_viz.read_and_modificate_input()
        comparison_path_template_2 = "{}/deseq_comp_{}_vs_{}_with_annotation_and_countings.csv".format(
            self._paths.deseq_extended_folder, conditions[1], conditions[0])
        comparison_2 = "{}_vs_{}".format(conditions[1], conditions[0])
        deseq_viz = DESeqViz(
            comparison_path_template_2, self._paths.viz_deseq_base_folder,
            self._args.padj_cutoff, comparison_2, self._args.alpha,
            self._args.color_sig, self._args.color_non_sig, self._args.shape,
            self._args.glyph_size)
        deseq_viz.read_and_modificate_input()
        

    
