import argparse
from reademptionlib.paths import Paths


def data_paths():
    global base_path
    global test_folder
    global test_files
    global test_lib_names
    global paths
    global folder_names
    global static_files
    parser = argparse.ArgumentParser()
    parser.add_argument("project_path", default="/tmp/test", nargs="?")
    args = parser.parse_args()
    args.project_path = "/tmp/test"
    paths = Paths(args)
    test_folder = "/tmp/test"
    test_files = ["foo.fa", "bar.fa"]
    test_lib_names = ["foo", "bar"]
    folder_names = [
        paths.input_folder,
        paths.output_folder,
        paths.align_report_folder,
        paths.raw_stat_data_folder,
        paths.read_fasta_folder,
        paths.ref_seq_folder,
        paths.annotation_folder,
        paths.read_alignment_index_folder,
        paths.read_alignments_folder,
        paths.read_processing_base_folder,
        paths.unaligned_reads_folder,
        paths.coverage_raw_folder,
        paths.coverage_tnoar_min_norm_folder,
        paths.coverage_tnoar_mil_norm_folder,
        paths.gene_quanti_base_folder,
        paths.gene_wise_quanti_combined_path]
    static_files = [
        paths.read_processing_stats_path,
        paths.read_alignments_stats_path,
        paths.read_file_stats,
        paths.read_alignment_stats_table_path,
        paths.ref_seq_file_stats,
        paths.index_path]
