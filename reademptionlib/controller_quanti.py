import concurrent.futures
import sys
import json
from reademptionlib.genewisequanti import GeneWiseOverview
from reademptionlib.genewisequanti import GeneWiseQuantification
from reademptionlib.helpers import Helpers
from reademptionlib.paths import Paths
from reademptionlib.rawstatdata import RawStatDataReader
from reademptionlib.visgenequanti import GeneQuantiVis


class QuantiController(object):

    def __init__(self, args):
        self._args = args
        self._paths = Paths(args)
        self._helpers = Helpers(args)

    def quantify_gene_wise(self):
        """Manage the counting of aligned reads per gene."""
        self._helpers.test_folder_existance(
            self._paths.required_gene_quanti_folders())
        norm_by_alignment_freq = True
        norm_by_overlap_freq = True
        if self._args.no_count_split_by_alignment_no:
            norm_by_alignment_freq = False
        if self._args.no_count_splitting_by_gene_no:
            norm_by_overlap_freq = False
        raw_stat_data_reader = RawStatDataReader()
        alignment_stats = [raw_stat_data_reader.read(
            self._paths.read_alignments_stats_path)]
        lib_names = sorted(list(alignment_stats[0].keys()))
        annotation_files = self._paths.get_annotation_files()
        self._paths.set_annotation_paths(annotation_files)
        was_paired_end_alignment = self._helpers.was_paired_end_alignment(
            lib_names)
        if not was_paired_end_alignment:
            self._paths.set_read_files_dep_file_lists_single_end(
                self._paths.get_read_files(), lib_names)
        else:
            self._paths.set_read_files_dep_file_lists_paired_end(
                self._paths.get_read_files(), lib_names)
        jobs = []
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=self._args.processes) as executor:
            for lib_name, read_alignment_path in zip(
                    lib_names, self._paths.read_alignment_bam_paths):
                jobs.append(executor.submit(
                    self._quantify_gene_wise, lib_name,
                    read_alignment_path, norm_by_alignment_freq,
                    norm_by_overlap_freq, annotation_files))
        # Evaluate thread outcome
        self._helpers.check_job_completeness(jobs)
        self._gene_quanti_create_overview(
            annotation_files, self._paths.annotation_paths, lib_names)
        self._vis_gene_quanti()

    def _quantify_gene_wise(
            self, lib_name, read_alignment_path,
            norm_by_alignment_freq,  norm_by_overlap_freq,
            annotation_files):
        """Perform the gene wise quantification for a given library."""
        gene_quanti_paths = [
            self._paths.gene_quanti_path(lib_name, annotation_file)
            for annotation_file in annotation_files]
        # Check if all output files for this library exist - if so
        # skip their creation
        if not any([self._helpers.file_needs_to_be_created(
                gene_quanti_path, quiet=True)
                for gene_quanti_path in gene_quanti_paths]):
            sys.stderr.write(
                "The file(s) %s exist(s). Skipping their/its generation.\n" %
                ", " .join(gene_quanti_paths))
            return
        gene_wise_quantification = GeneWiseQuantification(
            min_overlap=self._args.min_overlap,
            read_region=self._args.read_region,
            clip_length=self._args.clip_length,
            norm_by_alignment_freq=norm_by_alignment_freq,
            norm_by_overlap_freq=norm_by_overlap_freq,
            allowed_features_str=self._args.allowed_features,
            skip_antisense=self._args.skip_antisense,
            unique_only=self._args.unique_only)
        gene_wise_quantification.calc_overlaps_per_alignment(
            read_alignment_path, self._paths.annotation_paths)
        for annotation_file, annotation_path in zip(
                annotation_files, self._paths.annotation_paths):
            gene_wise_quantification.quantify(
                read_alignment_path, annotation_path,
                self._paths.gene_quanti_path(
                    lib_name, annotation_file), self._args.pseudocounts)

    def _gene_quanti_create_overview(
            self, annotation_files, annotation_paths, lib_names):
        """Create an overview table of all gene quantification for all libs."""
        strand_specific = True
        if self._args.non_strand_specific:
                    strand_specific = False
        gene_wise_overview = GeneWiseOverview(
            allowed_features_str=self._args.allowed_features,
            skip_antisense=self._args.skip_antisense,
            strand_specific=strand_specific)
        path_and_name_combos = {}
        for annotation_file, annotation_path in zip(
                annotation_files, annotation_paths):
            path_and_name_combos[annotation_path] = []
            for read_file in lib_names:
                path_and_name_combos[annotation_path].append(
                    [read_file, self._paths.gene_quanti_path(
                        read_file, annotation_file)])
        if self._helpers.file_needs_to_be_created(
                self._paths.gene_wise_quanti_combined_path):
            gene_wise_overview.create_overview_raw_countings(
                path_and_name_combos, lib_names,
                self._paths.gene_wise_quanti_combined_path)
        if self._helpers.file_needs_to_be_created(
                self._paths.gene_wise_quanti_combined_rpkm_path):
            gene_wise_overview.create_overview_rpkm(
                path_and_name_combos, lib_names,
                self._paths.gene_wise_quanti_combined_rpkm_path,
                self._libs_and_total_num_of_aligned_reads())
        if self._helpers.file_needs_to_be_created(
                self._paths.gene_wise_quanti_combined_tnoar_path):
            gene_wise_overview.create_overview_norm_by_tnoar(
                path_and_name_combos, lib_names,
                self._paths.gene_wise_quanti_combined_tnoar_path,
                self._libs_and_total_num_of_aligned_reads())

    def _libs_and_total_num_of_aligned_reads(self):
        """Read the total number of reads per library."""
        with open(self._paths
                  .read_alignments_stats_path) as read_aligner_stats_fh:
            read_aligner_stats = json.loads(read_aligner_stats_fh.read())
        return dict([(lib, values["stats_total"]["no_of_aligned_reads"])
                     for lib, values in read_aligner_stats.items()])

    def _vis_gene_quanti(self):
        """Generate plots based on the gene-wise read countings"""
        gene_quanti_vis = GeneQuantiVis(
            self._paths.gene_wise_quanti_combined_path,
            self._paths.get_lib_names_single_end() if not self._args.paired_end
            else self._paths.get_lib_names_paired_end(),
            self._paths.vis_gene_quanti_base_folder)
        gene_quanti_vis.parse_input_table()
        
