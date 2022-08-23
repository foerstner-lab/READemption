import sys
from collections import defaultdict
from collections import Counter
from functools import reduce
from reademptionlib.fasta import FastaParser
import pysam
from datetime import datetime

class ReadAlignerStats(object):
    def __init__(self, references_by_speies, paired_end=False, fragments=False):
        self.references_by_species = references_by_speies
        self.paired_end = paired_end
        self.fragments = fragments
        self.fasta_parser = FastaParser()
        if fragments:
            self.reads_or_fragments = "fragments"
        else:
            self.reads_or_fragments = "reads"
        """
        import itertools
        species = ["human", "virus", "bacteria"]
        cross_combos = []
        for combo_length in range(2, len(species) + 1):
            combos_of_length = list(itertools.combinations(species, combo_length))
            for combo_of_length in combos_of_length:
                cross_combos.append(combo_of_length)
        print(cross_combos)
        
        """

    def count(self, read_alignment_result_bam_path, unaligned_reads_path):
        self._stats = {}
        # Set up total stats
        self._stats["stats_total"] = defaultdict(float)
        self._init_stats_total()
        # Set up species stats
        self._stats["species_stats"] = defaultdict()
        for sp in self.references_by_species.keys():
            self._stats["species_stats"][sp] = defaultdict(float)
            self._init_species_dict(sp)
        print(f"readalignerstats count_aligned_reads_and_alignments start {datetime.now()}")
        self._count_aligned_reads_and_alignments(read_alignment_result_bam_path)
        print(f"readalignerstats count_aligned_reads_and_alignments stop {datetime.now()}")
        print(f"readalignerstats count_unaligned_reads start {datetime.now()}")
        self._count_unaligned_reads(unaligned_reads_path)
        print(f"readalignerstats count_unaligned_reads stop {datetime.now()}")
        return self._stats

    def _count_unaligned_reads(self, unaligned_read_paths):
        with open(unaligned_read_paths) as fasta_fh:
            self._stats["stats_total"][
                "no_of_unaligned_reads"
            ] = self._count_fasta_entries(fasta_fh)

    def _count_fasta_entries(self, fasta_fh):
        return reduce(
            lambda x, y: x + 1, self.fasta_parser.entries(fasta_fh), 0
        )

    def _count_aligned_reads_and_alignments(
        self, read_alignment_result_bam_path
    ):
        bam = pysam.Samfile(read_alignment_result_bam_path)
        bamfile = pysam.AlignmentFile(read_alignment_result_bam_path, "rb")
        indexed_bam = pysam.IndexedReads(bamfile)
        # build index in memory
        indexed_bam.build()

        stats_per_ref = defaultdict(dict)
        for ref_id in bam.references:
            # Set up reference stats
            self._init_counting_dict(stats_per_ref, ref_id)
        for entry in bam.fetch(until_eof=True):
            ref_id = bam.getrname(entry.tid)
            # Don't count the alignment if it is supplementary, to avoid
            # counting the same alignment multiple times
            if entry.is_supplementary:
                continue
            else:
                try:
                    self._count_alignment(
                        entry,
                        ref_id,
                        stats_per_ref,
                        indexed_bam,
                    )
                except KeyError:
                    sys.stderr.write(
                        "SAM entry with unspecified reference found! Stoping\n"
                    )
                    sys.exit(2)
        self._stats["stats_per_reference"] = stats_per_ref

    def _bam_to_sam(self, bam_path, sam_path):
        pysam.view("-ho{}".format(sam_path), bam_path, catch_stdout=False)

    def _init_stats_total(self):
        self._stats["stats_total"]["no_of_alignments"]
        self._stats["stats_total"][f"no_of_aligned_{self.reads_or_fragments}"]
        # self._stats["stats_total"]["fractions_of_aligned_reads"] # deprecated
        self._stats["stats_total"]["no_of_split_alignments"]
        self._stats["stats_total"][f"no_of_split_aligned_{self.reads_or_fragments}"]
        self._stats["stats_total"][f"no_of_uniquely_aligned_{self.reads_or_fragments}"]
        self._stats["stats_total"][f"no_of_multiple_aligned_{self.reads_or_fragments}"]
        self._stats["stats_total"]["alignment_length_and_freqs"] = defaultdict(
            int
        )
        self._stats["stats_total"][
            "no_of_hits_per_read_and_freqs"
        ] = defaultdict(int)
        self._stats["stats_total"][f"no_of_cross_aligned_{self.reads_or_fragments}"]

    def _init_species_dict(self, species):
        sp = species
        self._stats["species_stats"][sp]["no_of_alignments"]
        self._stats["species_stats"][sp][f"no_of_aligned_{self.reads_or_fragments}"]
        # self._stats["species_stats"][sp]["fractions_of_aligned_reads"] # deprecated
        self._stats["species_stats"][sp]["no_of_split_alignments"]
        self._stats["species_stats"][sp][f"no_of_split_aligned_{self.reads_or_fragments}"]
        self._stats["species_stats"][sp][f"no_of_uniquely_aligned_{self.reads_or_fragments}"]
        self._stats["species_stats"][sp][f"no_of_multiple_aligned_{self.reads_or_fragments}"]
        self._stats["species_stats"][sp][
            "alignment_length_and_freqs"
        ] = defaultdict(int)
        self._stats["species_stats"][sp][
            "no_of_hits_per_read_and_freqs"
        ] = defaultdict(int)
        self._stats["species_stats"][sp][f"no_of_cross_aligned_{self.reads_or_fragments}"]

    def _init_counting_dict(self, stats_per_ref, ref_id):
        stats_per_ref[ref_id] = defaultdict(float)
        stats_per_ref[ref_id]["no_of_alignments"]
        stats_per_ref[ref_id][f"no_of_aligned_{self.reads_or_fragments}"]
        # stats_per_ref[ref_id]["fractions_of_aligned_reads"] # deprecated
        stats_per_ref[ref_id]["no_of_split_alignments"]
        stats_per_ref[ref_id][f"no_of_split_aligned_{self.reads_or_fragments}"]
        stats_per_ref[ref_id][f"no_of_uniquely_aligned_{self.reads_or_fragments}"]
        stats_per_ref[ref_id][f"no_of_multiple_aligned_{self.reads_or_fragments}"]
        stats_per_ref[ref_id]["alignment_length_and_freqs"] = defaultdict(int)
        stats_per_ref[ref_id]["no_of_hits_per_read_and_freqs"] = defaultdict(
            int
        )

    def _count_alignment(
        self,
        entry,
        ref_id,
        stats_per_ref,
        indexed_bam,
    ):
        entry_tags_dict = dict(entry.tags)
        no_of_hits = entry_tags_dict["NH"]
        # check if alignment is unique or multiple
        if no_of_hits == 1:
            unique_alignment = True
        else:
            unique_alignment = False
        # check if alignment is split
        if "XH" in entry_tags_dict:
            split_alignment = True
        else:
            split_alignment = False
        # Add the alignment length frequencies to the chromosome
        stats_per_ref[ref_id]["alignment_length_and_freqs"][
            entry.reference_length
        ] += 1
        # Add to total stats if the alignment is primary (=not secondary)
        if not entry.is_secondary:
            # is primary alignment
            # get the species by ref_id
            for sp, ref_ids in self.references_by_species.items():
                if ref_id in ref_ids:
                    ref_sp = sp

            # Count number of hits and frequencies for the library
            self._stats["stats_total"]["no_of_hits_per_read_and_freqs"][
                no_of_hits
            ] += 1
            if unique_alignment and split_alignment:
                # Add to chromosome
                # Add to number of split aligned reads of chromosome
                stats_per_ref[ref_id][f"no_of_split_aligned_{self.reads_or_fragments}"] += 1
                # Add to number of aligned reads of chromosome
                stats_per_ref[ref_id][f"no_of_aligned_{self.reads_or_fragments}"] += 1
                # Add to number of alignments of chromosome
                stats_per_ref[ref_id]["no_of_alignments"] += 1
                # Count number of hits and frequencies of chromosome
                stats_per_ref[ref_id]["no_of_hits_per_read_and_freqs"][
                    no_of_hits
                ] += 1
                # Add to library
                # Add to number of split aligned reads of library
                self._stats["stats_total"][f"no_of_split_aligned_{self.reads_or_fragments}"] += 1
                # Add to number of aligned reads of library
                self._stats["stats_total"][f"no_of_aligned_{self.reads_or_fragments}"] += 1
                # Add to number of alignments of library
                self._stats["stats_total"]["no_of_alignments"] += 1
                # Add to species stats
                # Add to number of split aligned reads of species
                self._stats["species_stats"][ref_sp][
                    f"no_of_split_aligned_{self.reads_or_fragments}"
                ] += 1
                # Add to number of aligned reads of species
                self._stats["species_stats"][ref_sp][f"no_of_aligned_{self.reads_or_fragments}"] += 1
                # Add to number of alignments of species
                self._stats["species_stats"][ref_sp]["no_of_alignments"] += 1
                # Count number of hits and frequencies for the species
                self._stats["species_stats"][ref_sp][
                    "no_of_hits_per_read_and_freqs"
                ][no_of_hits] += 1

            # Count uniquely aligned read
            elif unique_alignment and not split_alignment:
                # Add to chromosome
                # Add to number of uniquely aligned reads of chromosome
                stats_per_ref[ref_id][f"no_of_uniquely_aligned_{self.reads_or_fragments}"] += 1
                # Add to number of aligned reads of chromosome
                stats_per_ref[ref_id][f"no_of_aligned_{self.reads_or_fragments}"] += 1
                # Add to number of alignments of chromosome
                stats_per_ref[ref_id]["no_of_alignments"] += 1
                # Count number of hits and frequencies of chromosome
                stats_per_ref[ref_id]["no_of_hits_per_read_and_freqs"][
                    no_of_hits
                ] += 1
                # Add to library
                # Add to number of uniquely aligned reads of library
                self._stats["stats_total"][f"no_of_uniquely_aligned_{self.reads_or_fragments}"] += 1
                # Add to number of aligned reads of library
                self._stats["stats_total"][f"no_of_aligned_{self.reads_or_fragments}"] += 1
                # Add to number of alignments of library
                self._stats["stats_total"]["no_of_alignments"] += 1
                # Add to species stats
                # Add to number of uniquely aligned reads of species
                self._stats["species_stats"][ref_sp][
                    f"no_of_uniquely_aligned_{self.reads_or_fragments}"
                ] += 1
                # Add to number of aligned reads of species
                self._stats["species_stats"][ref_sp][f"no_of_aligned_{self.reads_or_fragments}"] += 1
                # Add to number of alignments of species
                self._stats["species_stats"][ref_sp]["no_of_alignments"] += 1
                # Count number of hits and frequencies for the species
                self._stats["species_stats"][ref_sp][
                    "no_of_hits_per_read_and_freqs"
                ][no_of_hits] += 1

            # Count multiple aligned read
            elif not unique_alignment and not split_alignment:
                # Add to number of aligned reads of library
                self._stats["stats_total"][f"no_of_aligned_{self.reads_or_fragments}"] += 1
                # retrieve all alignments of the query
                alignments = indexed_bam.find(entry.query_name)
                # for paired end reads: select only the alignments of the current mate.
                if self.paired_end:
                    alignments = [alignment for alignment in alignments if alignment.is_read1 == entry.is_read1]
                # collect all reference names of alignments of query
                alignments_ref_seqs = []
                for alignment in alignments:
                    # do not collect supplementary alignments
                    if alignment.is_supplementary:
                        continue
                    alignments_ref_seqs.append(alignment.reference_name)
                for ref in alignments_ref_seqs:
                    # Add to number of alignments of chromosome
                    stats_per_ref[ref]["no_of_alignments"] += 1
                    # Add to number of alignments of library
                    self._stats["stats_total"]["no_of_alignments"] += 1

                # check if cross species aligned
                aligned_species = self._get_aligned_species(
                    alignments_ref_seqs, self.references_by_species
                )
                if len(set(aligned_species)) > 1:
                    # cross aligned
                    for ref_sp in set(aligned_species):
                        # Count number of hits and frequencies for the species
                        self._stats["species_stats"][ref_sp][
                            "no_of_hits_per_read_and_freqs"
                        ][no_of_hits] += 1
                        # Add to number of aligned reads to species
                        self._stats["species_stats"][ref_sp][
                            f"no_of_aligned_{self.reads_or_fragments}"
                        ] += 1
                        # Add to number of cross aligned reads to species
                        self._stats["species_stats"][ref_sp][
                            f"no_of_cross_aligned_{self.reads_or_fragments}"
                        ] += 1
                    for ref_sp in aligned_species:
                        # Add to number of alignments to species
                        self._stats["species_stats"][ref_sp][
                            "no_of_alignments"
                        ] += 1
                    # Add to number of cross aligned reads of library
                    self._stats["stats_total"][f"no_of_cross_aligned_{self.reads_or_fragments}"] += 1
                    for ref in set(alignments_ref_seqs):
                        # Add to number of cross aligned reads of chromosome.
                        # A set is used to ensure that a read that maps multiple
                        # times to the same chromosome is counted only once for
                        # each chromosome.
                        # Add to number of cross aligned reads of chromosome
                        stats_per_ref[ref][f"no_of_cross_aligned_{self.reads_or_fragments}"] += 1
                        # Add to number of aligned reads of chromosome
                        stats_per_ref[ref][f"no_of_aligned_{self.reads_or_fragments}"] += 1
                        # Count number of hits and frequencies of chromosome
                        stats_per_ref[ref]["no_of_hits_per_read_and_freqs"][
                            no_of_hits
                        ] += 1

                else:
                    # multiple aligned
                    ref_sp = aligned_species[0]
                    # Add to number of aligned reads of species
                    self._stats["species_stats"][ref_sp][
                        f"no_of_aligned_{self.reads_or_fragments}"
                    ] += 1
                    # Add to number of multiple aligned reads of species
                    self._stats["species_stats"][ref_sp][
                        f"no_of_multiple_aligned_{self.reads_or_fragments}"
                    ] += 1
                    # Add to number of multiple aligned reads of library
                    self._stats["stats_total"][
                        f"no_of_multiple_aligned_{self.reads_or_fragments}"
                    ] += 1

                    for ref in set(alignments_ref_seqs):
                        # Add to number of multiple aligned reads of chromosome.
                        # A set is used to ensure that a read that maps multiple
                        # times to the same chromosome is counted only once for
                        # each chromosome
                        stats_per_ref[ref][f"no_of_multiple_aligned_{self.reads_or_fragments}"] += 1
                        # Add to number of aligned reads of chromosome
                        stats_per_ref[ref][f"no_of_aligned_{self.reads_or_fragments}"] += 1
                        # Count number of hits and frequencies of chromosome
                        stats_per_ref[ref]["no_of_hits_per_read_and_freqs"][
                            no_of_hits
                        ] += 1
                    for ref_sp in aligned_species:
                        # Add to number of alignments to species
                        self._stats["species_stats"][ref_sp][
                            "no_of_alignments"
                        ] += 1
                    for ref_sp in set(aligned_species):
                        # Count number of hits and frequencies for the species
                        self._stats["species_stats"][ref_sp][
                            "no_of_hits_per_read_and_freqs"
                        ][no_of_hits] += 1

    def _get_aligned_species(self, alignment_ref_seqs, references_by_species):
        aligned_species = []
        for species, references in references_by_species.items():
            for ref_seq in alignment_ref_seqs:
                if ref_seq in references:
                    aligned_species.append(species)
        return aligned_species
