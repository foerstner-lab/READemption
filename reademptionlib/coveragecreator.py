import os
import sys
from reademptionlib.coveragecalculator import CoverageCalculator
from reademptionlib.wiggle import WiggleWriter


class CoverageCreator(object):
    def __init__(
        self,
        args,
        strands,
        coverage_files,
        species_references,
        count_cross_aligned_reads=False,
        cross_mapped_reads=None,
    ):
        self._args = args
        self.strands = strands
        self.coverage_files = coverage_files
        self._species_references = species_references
        self._count_cross_aligned_reads = count_cross_aligned_reads
        self._cross_mapped_reads = cross_mapped_reads

    def create_coverage_files_for_lib(
        self, lib_name, bam_path, no_of_aligned_reads, min_no_of_aligned_reads
    ):
        """Perform the coverage calculation for a given library."""
        # Don't create a coverage file, if no reads were aligned
        if no_of_aligned_reads == 0:
            return
        if self._all_coverage_file_exist():
            return
        read_count_splitting = True
        if self._args.skip_read_count_splitting:
            read_count_splitting = False
        coverage_calculator = CoverageCalculator(
            species_references=self._species_references,
            read_count_splitting=read_count_splitting,
            uniquely_aligned_only=self._args.unique_only,
            coverage_style=self._args.coverage_style,
            clip_length=self._args.clip_length,
            non_strand_specific=self._args.non_strand_specific,
            count_cross_aligned_reads=self._count_cross_aligned_reads,
            cross_mapped_reads=self._cross_mapped_reads,
        )
        (
            coverage_writers_raw,
            coverage_writers_tnoar_min_norm,
            coverage_writers_tnoar_mil_norm,
        ) = self._wiggle_writers(lib_name)
        for (
            ref_seq,
            coverages,
        ) in coverage_calculator.ref_seq_and_coverages(
            bam_path
        ):
            for strand in self.strands:
                coverage_writers_raw[strand].write_replicons_coverages(
                    ref_seq, coverages[strand]
                )
                coverage_writers_tnoar_min_norm[
                    strand
                ].write_replicons_coverages(
                    ref_seq,
                    coverages[strand],
                    factor=min_no_of_aligned_reads / no_of_aligned_reads,
                )
                coverage_writers_tnoar_mil_norm[
                    strand
                ].write_replicons_coverages(
                    ref_seq,
                    coverages[strand],
                    factor=1000000 / no_of_aligned_reads,
                )
        for strand in self.strands:
            coverage_writers_raw[strand].close_file()

    def _wiggle_writers(self, lib_name):
        """Write the calculated coverages to wiggle files."""
        coverage_writers_raw = dict(
            [
                (
                    strand,
                    WiggleWriter(
                        "%s_%s" % (lib_name, strand),
                        open(
                            self.coverage_files["wiggle_file_raw_path"][strand],
                            "w",
                        ),
                    ),
                )
                for strand in self.strands
            ]
        )
        coverage_writers_tnoar_min_norm = dict(
            [
                (
                    strand,
                    WiggleWriter(
                        "%s_%s" % (lib_name, strand),
                        open(
                            self.coverage_files[
                                "wiggle_file_tnoar_norm_min_path"
                            ][strand],
                            "w",
                        ),
                    ),
                )
                for strand in self.strands
            ]
        )
        coverage_writers_tnoar_mil_norm = dict(
            [
                (
                    strand,
                    WiggleWriter(
                        "%s_%s" % (lib_name, strand),
                        open(
                            self.coverage_files[
                                "wiggle_file_tnoar_norm_mil_path"
                            ][strand],
                            "w",
                        ),
                    ),
                )
                for strand in self.strands
            ]
        )
        return (
            coverage_writers_raw,
            coverage_writers_tnoar_min_norm,
            coverage_writers_tnoar_mil_norm,
        )

    def _all_coverage_file_exist(self):
        """Test the existance of all coverage files of a library"""
        files = []
        for strand in self.strands:
            files.append(self.coverage_files["wiggle_file_raw_path"][strand])
            files.append(
                self.coverage_files["wiggle_file_tnoar_norm_min_path"][strand]
            )

            files.append(
                self.coverage_files["wiggle_file_tnoar_norm_mil_path"][strand]
            )

        if not any(
            [self._file_needs_to_be_created(file, quiet=True) for file in files]
        ):
            sys.stderr.write(
                "The files %s exists. Skipping their generation.\n"
                % ", ".join(files)
            )
            return True
        return False

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
