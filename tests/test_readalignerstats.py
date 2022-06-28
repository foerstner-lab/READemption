import sys
import os
import unittest

sys.path.append(".")

from reademptionlib.readalignerstats import ReadAlignerStats
from collections import defaultdict


class ReadalignerStats(unittest.TestCase):
    maxDiff = None

    def test_count(self):
        library_one_aligned_path = (
            "tests/test_files/library_one_alignments_final.bam"
        )
        library_one_unaligned_path = "tests/test_files/library_one_unaligned.fa"
        references_by_species = {
            "human": [
                "GL000008.2",
                "chr1",
                "chr2",
                "chr3",
                "chr19",
                "chr20",
                "chr21",
                "chr22",
                "chr23",
                "chrX",
                "KQ031388.1",
            ],
            "staphylococcus": ["NC_007795.1"],
            "influenza": ["NC_007373.1", "NC_007372.1"],
        }
        read_aligner_stats = ReadAlignerStats(references_by_species)
        stats = read_aligner_stats.count(
            library_one_aligned_path, library_one_unaligned_path
        )
        expected_stats = {
            "species_stats": defaultdict(
                None,
                {
                    "human": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(int, {}),
                            "no_of_aligned_reads": 9.0,
                            "no_of_alignments": 11.0,
                            "no_of_cross_aligned_reads": 2.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {1: 6, 2: 1, 3: 2}
                            ),
                            "no_of_multiple_aligned_reads": 1.0,
                            "no_of_split_aligned_reads": 3.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 3.0,
                        },
                    ),
                    "influenza": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(int, {}),
                            "no_of_aligned_reads": 3.0,
                            "no_of_alignments": 3.0,
                            "no_of_cross_aligned_reads": 2.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {1: 1, 2: 1, 3: 1}
                            ),
                            "no_of_multiple_aligned_reads": 0.0,
                            "no_of_split_aligned_reads": 0.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 1.0,
                        },
                    ),
                    "staphylococcus": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(int, {}),
                            "no_of_aligned_reads": 5.0,
                            "no_of_alignments": 5.0,
                            "no_of_cross_aligned_reads": 3.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {1: 2, 2: 2, 3: 1}
                            ),
                            "no_of_multiple_aligned_reads": 0.0,
                            "no_of_split_aligned_reads": 0.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 2.0,
                        },
                    ),
                },
            ),
            "stats_per_reference": defaultdict(
                dict,
                {
                    "GL000008.2": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(int, {}),
                            "no_of_aligned_reads": 0.0,
                            "no_of_alignments": 0.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {}
                            ),
                            "no_of_multiple_aligned_reads": 0.0,
                            "no_of_split_aligned_reads": 0.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 0.0,
                        },
                    ),
                    "KQ031388.1": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(int, {}),
                            "no_of_aligned_reads": 0.0,
                            "no_of_alignments": 0.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {}
                            ),
                            "no_of_multiple_aligned_reads": 0.0,
                            "no_of_split_aligned_reads": 0.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 0.0,
                        },
                    ),
                    "NC_007372.1": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(
                                int, {80: 1}
                            ),
                            "no_of_aligned_reads": 1.0,
                            "no_of_alignments": 1.0,
                            "no_of_cross_aligned_reads": 1.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {3: 1}
                            ),
                            "no_of_multiple_aligned_reads": 0.0,
                            "no_of_split_aligned_reads": 0.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 0.0,
                        },
                    ),
                    "NC_007373.1": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(
                                int, {80: 2}
                            ),
                            "no_of_aligned_reads": 2.0,
                            "no_of_alignments": 2.0,
                            "no_of_cross_aligned_reads": 1.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {2: 1, 1: 1}
                            ),
                            "no_of_multiple_aligned_reads": 0.0,
                            "no_of_split_aligned_reads": 0.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 1.0,
                        },
                    ),
                    "NC_007795.1": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(
                                int, {80: 5}
                            ),
                            "no_of_aligned_reads": 5.0,
                            "no_of_alignments": 5.0,
                            "no_of_cross_aligned_reads": 3.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {1: 2, 2: 2, 3: 1}
                            ),
                            "no_of_multiple_aligned_reads": 0.0,
                            "no_of_split_aligned_reads": 0.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 2.0,
                        },
                    ),
                    "chr1": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(
                                int, {60: 2, 80: 1}
                            ),
                            "no_of_aligned_reads": 3.0,
                            "no_of_alignments": 3.0,
                            "no_of_cross_aligned_reads": 1.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {1: 2, 3: 1}
                            ),
                            "no_of_multiple_aligned_reads": 0.0,
                            "no_of_split_aligned_reads": 0.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 2.0,
                        },
                    ),
                    "chr19": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(
                                int, {39: 1, 651: 1}
                            ),
                            "no_of_aligned_reads": 2.0,
                            "no_of_alignments": 2.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {1: 2}
                            ),
                            "no_of_multiple_aligned_reads": 0.0,
                            "no_of_split_aligned_reads": 2.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 0.0,
                        },
                    ),
                    "chr2": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(
                                int, {40: 2, 60: 1}
                            ),
                            "no_of_aligned_reads": 2.0,
                            "no_of_alignments": 3.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {1: 1, 3: 1}
                            ),
                            "no_of_multiple_aligned_reads": 1.0,
                            "no_of_split_aligned_reads": 0.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 1.0,
                        },
                    ),
                    "chr20": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(
                                int, {80: 1}
                            ),
                            "no_of_aligned_reads": 1.0,
                            "no_of_alignments": 1.0,
                            "no_of_cross_aligned_reads": 1.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {2: 1}
                            ),
                            "no_of_multiple_aligned_reads": 0.0,
                            "no_of_split_aligned_reads": 0.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 0.0,
                        },
                    ),
                    "chr21": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(int, {}),
                            "no_of_aligned_reads": 0.0,
                            "no_of_alignments": 0.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {}
                            ),
                            "no_of_multiple_aligned_reads": 0.0,
                            "no_of_split_aligned_reads": 0.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 0.0,
                        },
                    ),
                    "chr22": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(
                                int, {120: 1}
                            ),
                            "no_of_aligned_reads": 1.0,
                            "no_of_alignments": 1.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {1: 1}
                            ),
                            "no_of_multiple_aligned_reads": 0.0,
                            "no_of_split_aligned_reads": 1.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 0.0,
                        },
                    ),
                    "chr23": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(int, {}),
                            "no_of_aligned_reads": 0.0,
                            "no_of_alignments": 0.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {}
                            ),
                            "no_of_multiple_aligned_reads": 0.0,
                            "no_of_split_aligned_reads": 0.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 0.0,
                        },
                    ),
                    "chr3": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(
                                int, {40: 1}
                            ),
                            "no_of_aligned_reads": 1.0,
                            "no_of_alignments": 1.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {3: 1}
                            ),
                            "no_of_multiple_aligned_reads": 1.0,
                            "no_of_split_aligned_reads": 0.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 0.0,
                        },
                    ),
                    "chrX": defaultdict(
                        float,
                        {
                            "alignment_length_and_freqs": defaultdict(int, {}),
                            "no_of_aligned_reads": 0.0,
                            "no_of_alignments": 0.0,
                            "no_of_hits_per_read_and_freqs": defaultdict(
                                int, {}
                            ),
                            "no_of_multiple_aligned_reads": 0.0,
                            "no_of_split_aligned_reads": 0.0,
                            "no_of_split_alignments": 0.0,
                            "no_of_uniquely_aligned_reads": 0.0,
                        },
                    ),
                },
            ),
            "stats_total": defaultdict(
                float,
                {
                    "alignment_length_and_freqs": defaultdict(int, {}),
                    "no_of_aligned_reads": 13.0,
                    "no_of_alignments": 19.0,
                    "no_of_cross_aligned_reads": 3.0,
                    "no_of_hits_per_read_and_freqs": defaultdict(
                        int, {1: 9, 2: 2, 3: 2}
                    ),
                    "no_of_multiple_aligned_reads": 1.0,
                    "no_of_split_aligned_reads": 3.0,
                    "no_of_split_alignments": 0.0,
                    "no_of_unaligned_reads": 1,
                    "no_of_uniquely_aligned_reads": 6.0,
                },
            ),
        }
        self.assertDictEqual(expected_stats, stats)
