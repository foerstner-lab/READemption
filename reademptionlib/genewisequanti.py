import csv
from reademptionlib.gff3 import Gff3Parser
import pysam
import pandas as pd


class GeneWiseQuantification(object):
    def __init__(
        self,
        references_by_species,
        min_overlap=1,
        read_region="global",
        clip_length=11,
        norm_by_alignment_freq=True,
        norm_by_overlap_freq=True,
        allowed_features_str=None,
        add_antisense=False,
        antisense_only=False,
        strand_specific=True,
        unique_only=False,
        count_cross_aligned_reads=False,
        crossmapped_reads=None,
    ):
        """
        - normalize_by_alignment: consider that some reads are aligned at
          more than one location and only count fractions
        - normalize_by_overlapping_genes: consider that some alignment
          overlap with more than on gene

        """
        self._references_by_species = references_by_species
        self._min_overlap = min_overlap
        self._read_region = read_region
        self._clip_length = clip_length
        self._norm_by_alignment_freq = norm_by_alignment_freq
        self._norm_by_overlap_freq = norm_by_overlap_freq
        self._allowed_features = _allowed_features(allowed_features_str)
        self._add_antisense = add_antisense
        self._antisense_only = antisense_only
        self._strand_specific = strand_specific
        self._unique_only = unique_only
        self._count_cross_aligned_reads = count_cross_aligned_reads
        self._crossmapped_reads = crossmapped_reads

    def calc_overlaps_per_alignment(
        self, read_alignment_path, annotation_paths
    ):
        """Calculate for each alignment the number of genes it
        overlaps. This has to be done globally i.e. for all annotation
        files combined in one dictionary.
        """
        gff3_parser = Gff3Parser()
        self.alignments_and_no_of_overlaps = {}
        for annotation_path in annotation_paths:
            annotation_name = annotation_path.split("/")[-1]
            sam = pysam.Samfile(read_alignment_path)
            for entry in gff3_parser.entries(
                open(annotation_path), annotation_name
            ):
                if _entry_to_use(entry, self._allowed_features) is False:
                    continue
                for alignment in self._overlapping_alignments(sam, entry):
                    alignment_id = self._alignment_id(alignment)
                    # check if species cross aligned reads should be counted
                    if not self._count_cross_aligned_reads:
                        # skip alignment if it is cross aligned
                        if alignment.qname in self._crossmapped_reads:
                            continue
                    self.alignments_and_no_of_overlaps.setdefault(
                        alignment_id, 0
                    )
                    self.alignments_and_no_of_overlaps[alignment_id] += 1

    def quantify(
        self,
        read_alignment_path,
        annotation_path,
        output_path,
        pseudocounts=False,
    ):
        self._quantify(
            read_alignment_path,
            annotation_path,
            output_path,
            self._fraction_calc_method(),
            pseudocounts,
        )

    def _quantify(
        self,
        read_alignment_path,
        annotation_path,
        output_path,
        fraction_calc_method,
        pseudocounts=False,
    ):

        sam = pysam.Samfile(read_alignment_path)
        gff3_parser = Gff3Parser()
        output_fh = open(output_path, "w")
        output_fh.write(
            "#"
            + "\t".join(_gff_field_descriptions() + ["sense", "antisense"])
            + "\n"
        )
        annotation_name = annotation_path.split("/")[-1]
        for entry in gff3_parser.entries(
            open(annotation_path), annotation_name
        ):
            if _entry_to_use(entry, self._allowed_features) is False:
                continue
            if pseudocounts is False:
                sum_sense = 0
                sum_antisense = 0
            else:
                sum_sense = 1
                sum_antisense = 1
            for alignment in self._overlapping_alignments(sam, entry):
                # check if species cross aligned reads should be counted
                if not self._count_cross_aligned_reads:
                    # skip alignment if it is cross aligned
                    if alignment.qname in self._crossmapped_reads:
                        continue
                fraction = fraction_calc_method(alignment)
                if self._same_strand(entry, alignment):
                    sum_sense += fraction
                else:
                    sum_antisense += fraction
            output_fh.write(
                str(entry)
                + "\t"
                + str(sum_sense)
                + "\t"
                + str(sum_antisense)
                + "\n"
            )

    def _same_strand(self, entry, alignment):
        assert entry.strand in ["+", "-"]
        if alignment.is_read2 is False:
            if (entry.strand == "+" and alignment.is_reverse is False) or (
                entry.strand == "-" and alignment.is_reverse is True
            ):
                return True
        # Mate pair for paired end sequencing
        elif alignment.is_read2 is True:
            if (entry.strand == "+" and alignment.is_reverse is True) or (
                entry.strand == "-" and alignment.is_reverse is False
            ):
                return True
        return False

    def _fraction_calc_method(self):
        if self._norm_by_alignment_freq and self._norm_by_overlap_freq:
            return self._fraction_norm_by_alignment_and_overlap
        elif self._norm_by_alignment_freq and not self._norm_by_overlap_freq:
            return self._fraction_norm_by_alignment
        elif not self._norm_by_alignment_freq and self._norm_by_overlap_freq:
            return self._fraction_norm_by_overlap
        return self._fraction_calc_constant_one

    def _alignment_tags(self, alignment):
        return dict(alignment.tags)

    def _fraction_calc_constant_one(self, alignment):
        return 1.0

    def _fraction_norm_by_alignment_and_overlap(self, alignment):
        alignment_tags = self._alignment_tags(alignment)
        return (
            1.0
            / float(
                self.alignments_and_no_of_overlaps[
                    self._alignment_id(alignment)
                ]
            )
            / float(alignment_tags["NH"])  # no. of alignments of read
        )

    def _fraction_norm_by_alignment(self, alignment):
        alignment_tags = self._alignment_tags(alignment)
        return 1.0 / float(alignment_tags["NH"])  # no. of alignments of read

    def _fraction_norm_by_overlap(self, alignment):
        alignment_tags = self._alignment_tags(alignment)
        return 1.0 / float(
            self.alignments_and_no_of_overlaps[self._alignment_id(alignment)]
        )

    def _overlapping_alignments(self, sam, entry):
        # The substraction of 1 from the start is necessary to perform
        # this correctly (checked in IGB, IGV and the unit testings).
        for alignment in sam.fetch(
            reference=entry.seq_id, start=entry.start - 1, end=entry.end
        ):
            # 1-based alignment coordinates
            start = alignment.pos + 1
            end = alignment.aend
            if self._read_region == "first_base_only":
                if (alignment.is_reverse is False) and (
                    (start < entry.start) or (start > entry.end)
                ):
                    continue
                if (alignment.is_reverse is True) and (
                    (end < entry.start) or (end > entry.end)
                ):
                    continue
            elif self._read_region == "last_base_only":
                if (alignment.is_reverse is False) and (
                    (end < entry.start) or (end > entry.end)
                ):
                    continue
                if (alignment.is_reverse is True) and (
                    (start < entry.start) or (start > entry.end)
                ):
                    continue
            elif self._read_region == "centered":
                if (
                    _get_overlap(
                        start + self._clip_length,
                        end - self._clip_length,
                        entry.start,
                        entry.end,
                    )
                    < self._min_overlap
                ):
                    continue
            else:
                if (
                    alignment.get_overlap(entry.start - 1, entry.end)
                    < self._min_overlap
                ):
                    continue
            if (
                not self._add_antisense
                and not self._antisense_only
                and self._strand_specific
            ):
                if not self._same_strand(entry, alignment):
                    continue
            if self._antisense_only:
                if self._same_strand(entry, alignment):
                    continue
            if self._unique_only:
                if dict(alignment.tags)["NH"] != 1:
                    continue
            yield (alignment)

    def _alignment_id(self, alignment):
        return ":".join(
            [
                str(alignment.tid),
                alignment.qname,
                str(alignment.flag),
                str(alignment.pos),
                str(alignment.aend),
            ]
        )

    def _values_to_gene_key(self, seq_id, feature, start, end, strand):
        return "|".join(
            [str(val) for val in [seq_id, feature, start, end, strand]]
        )


class GeneWiseOverview(object):
    def __init__(
        self,
        allowed_features_str=None,
        add_antisense=False,
        antisense_only=False,
        strand_specific=True,
    ):
        self._allowed_features = _allowed_features(allowed_features_str)
        self._add_antisense = add_antisense
        self._antisense_only = antisense_only
        self._strand_specific = strand_specific

    def create_overview_raw_countings(
        self, path_and_name_combos, read_files, overview_path
    ):
        self._create_overview(path_and_name_combos, read_files, overview_path)

    def create_overview_rpkm(
        self, path_and_name_combos, read_files, overview_path, libs_and_tnoar
    ):
        self._create_overview(
            path_and_name_combos,
            read_files,
            overview_path,
            normalization="RPKM",
            libs_and_tnoar=libs_and_tnoar,
        )

    def create_overview_norm_by_tnoar(
        self, path_and_name_combos, read_files, overview_path, libs_and_tnoar
    ):
        self._create_overview(
            path_and_name_combos,
            read_files,
            overview_path,
            normalization="TNOAR",
            libs_and_tnoar=libs_and_tnoar,
        )

    def create_overview_tpm(
        self, gene_wise_quanti_combined_path, gene_wise_quanti_combined_tpm_path
    ):
        gene_quanti = pd.read_csv(gene_wise_quanti_combined_path, sep="\t")
        # the libs are starting at column 11
        libs = gene_quanti.columns.to_list()[10:]
        gene_quanti_tpm = self._calculate_tpm(gene_quanti, libs)
        gene_quanti_tpm.to_csv(
            gene_wise_quanti_combined_tpm_path, sep="\t", index=False
        )

    def _calculate_tpm(self, gene_quanti, libs) -> pd.DataFrame:
        """
        :param gene_quanti: a pandas data frame generated from the gene wise quantification
        table containing the raw reads
        :param libs: a list of library names extracted from the gene wise quantification table
        :return: a pandas data frame containing TPM values instead of raw read counts

        Formula to calculate TPM (transcripts per million) from
        "Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples",
        Günter P. Wagner, Koryu Kin & Vincent J. Lynch,
        DOI: 10.1007/s12064-012-0162-3

                r_g x rl x 1000000
        TPM  = ────────────────────
                    fl_g x T
        where
          r_g = number of reads that map to a gene
          rl = read length i.e., the average number of nucleotides mapped per read
          fl_g = feature length or length of the gene
          T is the total number of transcripts sampled in a sequencing run and is calculated as follows:
                ___
                ╲     r_g x rl
          T =   ╱    ─────────
                ‾‾‾     fl_g
               g e G

        The Formula can be simplified (by excluding the read length rl) to:

                r_g x 1000000
        TPM  = ──────────────
                   fl_g x A
        where
                ___
                ╲     r_g
          A =   ╱    ────
                ‾‾‾   fl_g
               g e G
         The simplified formula is implemented below
        """
        for lib in libs:
            gene_quanti[lib] = gene_quanti[lib].astype(float)
            if (gene_quanti[lib] == 0).all():
                print(
                    f"Warning: Calculating TPM values for genes that have no "
                    f"other values than zero is not possible. Skipping the "
                    f"creation of the TPM gene quantification for library {lib}."
                )
                gene_quanti.drop(lib, inplace=True, axis=1)
                continue
            # calculate A
            gene_quanti["transcript_count"] = gene_quanti.apply(
                lambda df: (float(df[lib]))
                / (int(df["End"]) - int(df["Start"]) + 1),
                axis=1,
            )
            A = gene_quanti["transcript_count"].sum()
            # calculate TPM per gene and replace the raw read counts in the gene quanti table
            gene_quanti[lib] = gene_quanti.apply(
                lambda df: (float(df[lib]) * 1000000)
                / ((int(df["End"]) - int(df["Start"]) + 1) * A),
                axis=1,
            )
            gene_quanti.drop("transcript_count", inplace=True, axis=1)
        return gene_quanti

    def _create_overview(
        self,
        path_and_name_combos,
        read_files,
        overview_path,
        normalization=None,
        libs_and_tnoar=None,
    ):
        output_fh = open(overview_path, "w")
        # Write header
        output_fh.write(
            "\t".join(
                [
                    "Orientation of counted reads relative to the strand "
                    "location of the annotation"
                ]
                + _gff_field_descriptions()
                + read_files
            )
            + "\n"
        )
        if self._strand_specific and not self._antisense_only:
            self._add_to_overview(
                path_and_name_combos,
                "sense",
                9,
                output_fh,
                normalization,
                libs_and_tnoar,
            )
        if self._add_antisense or self._antisense_only:
            self._add_to_overview(
                path_and_name_combos,
                "anti-sense",
                10,
                output_fh,
                normalization,
                libs_and_tnoar,
            )
        if not self._strand_specific:
            self._add_to_overview_strand_unspecific(
                path_and_name_combos,
                "sense_and_antisense",
                9,
                10,
                output_fh,
                normalization,
                libs_and_tnoar,
            )

    def _add_to_overview(
        self,
        path_and_name_combos,
        direction,
        column,
        output_fh,
        normalization=None,
        libs_and_tnoar=None,
    ):
        gff3_parser = Gff3Parser()
        for annotation_path in sorted(path_and_name_combos.keys()):
            table_columns = []
            entries = []
            seq_lengths = []
            annotation_name = annotation_path.split("/")[-1]
            for entry in gff3_parser.entries(
                open(annotation_path), annotation_name
            ):
                if _entry_to_use(entry, self._allowed_features) is False:
                    continue
                entries.append(direction + "\t" + str(entry))
                seq_lengths.append(entry.end - entry.start + 1)
            table_columns.append(entries)
            for read_file, gene_quanti_path in path_and_name_combos[
                annotation_path
            ]:
                reader = csv.reader(open(gene_quanti_path), delimiter="\t")
                next(reader)  # skip first line
                if normalization == "RPKM":
                    table_columns.append(
                        [
                            self._rpkm(
                                row[column], length, libs_and_tnoar[read_file]
                            )
                            for row, length in zip(reader, seq_lengths)
                        ]
                    )
                elif normalization == "TNOAR":
                    table_columns.append(
                        [
                            self._norm_by_tnoar(
                                row[column], libs_and_tnoar[read_file]
                            )
                            for row, length in zip(reader, seq_lengths)
                        ]
                    )
                else:
                    table_columns.append([row[column] for row in reader])
            # Generate a table by rotating the column list
            table = zip(*table_columns)
            for row in table:
                output_fh.write("\t".join(row) + "\n")

    def _add_to_overview_strand_unspecific(
        self,
        path_and_name_combos,
        direction,
        column1,
        column2,
        output_fh,
        normalization=None,
        libs_and_tnoar=None,
    ):
        gff3_parser = Gff3Parser()
        for annotation_path in sorted(path_and_name_combos.keys()):
            table_columns = []
            entries = []
            seq_lengths = []
            annotation_name = annotation_path.split("/")[-1]
            for entry in gff3_parser.entries(
                open(annotation_path), annotation_name
            ):
                if _entry_to_use(entry, self._allowed_features) is False:
                    continue
                entries.append(direction + "\t" + str(entry))
                seq_lengths.append(entry.end - entry.start + 1)
            table_columns.append(entries)
            for read_file, gene_quanti_path in path_and_name_combos[
                annotation_path
            ]:
                reader = csv.reader(open(gene_quanti_path), delimiter="\t")
                next(reader)  # skip first line
                if normalization == "RPKM":
                    table_columns.append(
                        [
                            self._rpkm(
                                str(float(row[column1]) + float(row[column2])),
                                length,
                                libs_and_tnoar[read_file],
                            )
                            for row, length in zip(reader, seq_lengths)
                        ]
                    )
                elif normalization == "TNOAR":
                    table_columns.append(
                        [
                            self._norm_by_tnoar(
                                str(float(row[column1]) + float(row[column2])),
                                libs_and_tnoar[read_file],
                            )
                            for row, length in zip(reader, seq_lengths)
                        ]
                    )
                else:
                    table_columns.append(
                        [
                            str(float(row[column1]) + float(row[column2]))
                            for row in reader
                        ]
                    )
            # Generate a table by rotating the column list
            table = zip(*table_columns)
            for row in table:
                output_fh.write("\t".join(row) + "\n")

    def _rpkm(self, counting, length, total_no_of_aligned_reads):
        """
        Formula in Supplemenatary Material S1 of
        http://www.nature.com/nmeth/journal/v5/n7/full/nmeth.1226.html

        R = (10^9 * C) / (N * L)

        with C = is the number of mappable reads that fell onto the gene
             N = total number of mappable read
             L = length of the gene

        """
        return str(
            float(counting)
            * float(10**9)
            / (float(total_no_of_aligned_reads) * float(length))
        )

    def _norm_by_tnoar(self, counting, total_no_of_aligned_reads):
        return str(float(counting) / float(total_no_of_aligned_reads))


def _entry_to_use(entry, allowed_features):
    if allowed_features is None:
        return True
    if entry.feature in allowed_features:
        return True
    return False


def _allowed_features(allowed_features_str):
    if allowed_features_str is None:
        return None
    else:
        return [feature.strip() for feature in allowed_features_str.split(",")]


def _gff_field_descriptions():
    return [
        "Sequence name",
        "Source",
        "Feature",
        "Start",
        "End",
        "Score",
        "Strand",
        "Frame",
        "Attributes",
    ]


def _get_overlap(alignment_start, alignment_end, feature_start, feature_end):
    return max(
        0,
        min(alignment_end, feature_end)
        - max(alignment_start, feature_start)
        + 1,
    )
