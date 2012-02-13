#!/usr/bin/env python

# Copyright (c) 2011, Konrad Foerstner <konrad@foerstner.org>
#
# Permission to use, copy, modify, and/or distribute this software for
# any purpose with or without fee is hereby granted, provided that the
# above copyright notice and this permission notice appear in all
# copies.
#
# THE SOFTWARE IS PROVIDED 'AS IS' AND THE AUTHOR DISCLAIMS ALL
# WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
# AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL
# DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA
# OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
# TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.

"""Combines overlap countings of read mappings and annotations of different files.

usage: build_annotation_table_with_read_countings.py [-h] [-m MIN_OVERLAP]
                                                     [-d STRAND_ORIENTATION]
                                                     [-n NORMALIZATION_FACTORS]
                                                     [-r]
                                                     ANNOTATION_FILE
                                                     OVERLAP_FILE
                                                     [OVERLAP_FILE ...]

positional arguments:
  ANNOTATION_FILE       An annotation table files (NCBI style).
  OVERLAP_FILE          Overlap file(s) produced by
                        sam_hit_annotation_mapping.py.

optional arguments:
  -h, --help            show this help message and exit
  -m MIN_OVERLAP, --min_overlap MIN_OVERLAP
                        minimal overlap needed to count a match (in bp).
  -p MIN_READ_PERCENTAGE, --min_read_percentage MIN_READ_PERCENTAGE
                        the minimal percentage value of the read mapping
                        length, that has to overlap with an annotation to be
                        counted. E.g. if a read mapping is 80 bp long and
                        overlaps with a annotation in 20 bp the percentage
                        value is 25.
  -d STRAND_ORIENTATION, --strand_orientation STRAND_ORIENTATION
                        direction of the sequence in in comparison to the
                        annotation. There are three posibilities: s: sense; a:
                        antisense; b: both. Default is s.
  -n NORMALIZATION_FACTORS [NORMALIZATION_FACTORS ...], 
                        --normalization_factors NORMALIZATION_FACTORS [NORMALIZATION_FACTORS ...]
                        They are used to normalized the countings (e.g. by
                        total number of mapped reads).
  -r, --rpkm            calculate RPKM instead of raw countings. Requires -n
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        An output file. If not given the output is written to
                        the standard output.
"""
__description__ = ("Combines overlap countings of read mappings and "
                   "annotations of different files.")
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2011 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = "0.1"

import sys
from argparse import ArgumentParser

def main():
    """
    Run, forest, run!
    """
    python_version = sys.version.split()[0]
    if not python_version[0] == "3":
        sys.stdout.write('Error! Please us Python 3.x!\n')
        sys.exit(2)
    parser = create_arg_parser()
    args = parser.parse_args()
    annotation_mapping_table_builder = AnnotationMappingTableBuilder(
        args.annotation_file, args.overlap_files, min_overlap=args.min_overlap,
        min_read_percentage=args.min_read_percentage, 
        output_file=args.output_file, rpkm=args.rpkm,
        strand_orientation=args.strand_orientation, 
        normalization_factors=args.normalization_factors,
        count_nucleotides=args.count_nucleotides)
    annotation_mapping_table_builder.check_input()
    annotation_mapping_table_builder.read_annotation_mapping_files()
    annotation_mapping_table_builder.read_annotation_file_and_print_output()
    annotation_mapping_table_builder.sum_up_cells()

def create_arg_parser():
    parser = ArgumentParser(description=__description__)
    parser.add_argument(
        "annotation_file", metavar='ANNOTATION_FILE', 
        help="An annotation table files (NCBI style).")
    parser.add_argument(
        "overlap_files", metavar='OVERLAP_FILE', nargs="+", help="Overlap "
        "file(s) produced by sam_hit_annotation_mapping.py.")
    parser.add_argument(
        "-m", "--min_overlap", default=1, action="store", 
        dest="min_overlap", help="minimal overlap needed to count a match "
        "(in bp).")
    parser.add_argument(
        "-p", "--min_read_percentage", default=0, action="store", 
        dest="min_read_percentage", help="the minimal percentage value of "
        "the read mapping length, that has to overlap with an annotation "
        "to be counted. E.g. if a read mapping is 80 bp long and overlaps "
        "with a annotation in 20 bp the percentage value is 25.")
    parser.add_argument(
        "-d", "--strand_orientation", default="s", action="store", 
        dest="strand_orientation", help="direction of the sequence in in "
        "comparison to the annotation. There are three posibilities: "
        "s: sense; a: antisense; b: both. Default is s.")
    parser.add_argument(
        "-n", "--normalization_factors", default=None, action="store", 
        dest="normalization_factors", nargs="+", 
        help="They are used to normalized the countings (e.g. by total "
        "number of mapped reads).")
    parser.add_argument("-r", "--rpkm", default=False, action="store_true", 
                        dest="rpkm", help="calculate RPKM instead of raw " +
                        "countings. Requires -n")
    parser.add_argument("-b", "--count_nucleotides", default=False, 
                        action="store_true", dest="count_nucleotides", 
                        help="count mapped nucleotides instead of mapped " +
                        "reads.")
    parser.add_argument("-o", "--output_file", default=None, action="store", 
                        dest="output_file", help="An output file. If not " +
                        "given the output is written to the standard output.")
    return(parser)

class AnnotationMappingTableBuilder(object):
    """Combines overlap countings of read mappings and annotations of 
    different files.

    It reads output produced by sam_hit_annotation_mapping.py.
    """

    def __init__(self, annotation_file, overlap_files, min_overlap=1, 
                 min_read_percentage=0, output_file=None, 
                 strand_orientation="s", normalization_factors=None, 
                 rpkm=False, count_nucleotides=False):
        self.annotation_file = annotation_file
        self.annotation_mapping_files = overlap_files
        self.min_overlap = int(min_overlap)
        self.min_read_percentage = float(min_read_percentage)
        self.normalization_factors = normalization_factors
        self.count_nucleotides = count_nucleotides
        self.rpkm = rpkm
        if strand_orientation not in ["s", "a", "b"]:
            sys.stdout.write(
                "Error! Invalid value for the strand orientation option.")
            sys.exit(2)
        self.strand_orientation = strand_orientation
        if output_file:
            self.output_fh = open(output_file, "w")
        else:
            self.output_fh = sys.stdout

        # Create for every file a dictionary that contains the
        # number of read mapping annoations overlaps countings per
        # gene.
        self.mapping_files_and_annotation_counting = {}

    def check_input(self):
        """Test consistency of different parameters."""
        if self.normalization_factors:
            # There must be one normalization factor per read mapping
            # annotion overlap file.
            if len(self.annotation_mapping_files) != len(
                self.normalization_factors):
                sys.stdout.write(
                    "Error! Number of normalization factor is not "
                    "equal the number of annotation mapping files.")
                sys.stdout.write(
                    "Normalisation factors: %s" % self.normalization_factors)
                sys.stdout.write("Files: %s" % self.annotation_mapping_files)
                sys.exit(2)
        if self.rpkm and not self.normalization_factors:
            sys.stdout.write("If -r is set, -n must be set, too. See -h.\n")
            sys.exit(2)
        if self.rpkm and self.count_nucleotides:
            sys.stdout.write("Choose either RPKM or nucleotide counting. "
                             "Not both.\n")
            sys.exit(2)

    def read_annotation_mapping_files(self):
        """Read the annotation mapping files an perform the countings."""
        for mapping_file in self.annotation_mapping_files:
            self.mapping_files_and_annotation_counting[mapping_file] = {}
            self._read_annotation_mapping_file(mapping_file)
            
    def read_annotation_file_and_print_output(self):
        """Read the annotation file and generate the output."""
        self.output_fh.write(self._parameter_dump())
        self.output_fh.write(self._headline())
        mapping_files = sorted(self.mapping_files_and_annotation_counting.keys())
        for line in open(self.annotation_file):
            split_line = line.split("\t")
            if '..' not in split_line[0] or len(split_line) != 9:
                continue
            entry = self._parse_annotation_line(line)
            key = self._annotation_entry_key(entry)
            countings = [self.mapping_files_and_annotation_counting[
                        mapping_file].get(key, 0)
                         for mapping_file in mapping_files]
            if self.rpkm and self.normalization_factors:
                countings = self._rpkm_normalized_countings(
                    countings, entry)
            elif self.normalization_factors:
                countings = self._normalized_countings(countings)
            countings = [str(counting) for counting in countings]
            self.output_fh.write(
                line[:-1] + "\t" + "\t".join(countings) + "\n")

    def sum_up_cells(self):
        """Sum up all cells for each mapping file"""
        mapping_files_and_sum = {}
        for mapping_file, annoations_and_coutings in (
            self.mapping_files_and_annotation_counting.items()):
            mapping_files_and_sum[mapping_file.split("/")[-1]] = sum(
                annoations_and_coutings.values())
        return(mapping_files_and_sum)

    def _parameter_dump(self):
        """Retur a string with selected parameters to be printed."""
        min_read_percentage = "-"
        if not self.min_read_percentage:
            min_read_percentage = self.min_read_percentage
        strand_orientation = {
            "s": "sense", "a" : "antisense", "b" : "sense or antisense"
            }[self.strand_orientation]
        return("# Minimal overlap [bp]: %s\n" 
               "# Minimal read mapping length percentage (%%): %s\n" 
               "# Required strand orientation: %s\n\n" % (
                str(self.min_overlap), str(min_read_percentage), 
                strand_orientation))

    def _normalized_countings(self, read_countings):
        """Normalize the countings by using the normalization factors."""
        return(
            [str(float(counting)/float(normalization_factor)) 
             for counting, normalization_factor
             in  zip(read_countings, self.normalization_factors)])

    def _rpkm_normalized_countings(self, read_countings, entry):
        """Normalized the counting to RPKM values."""
        gene_length = self._gene_length(entry)
        return(
            [str(self._calc_rpkm(counting, normalization_factor, gene_length))
             for counting, normalization_factor 
             in zip(read_countings, self.normalization_factors)])

    def _headline(self):
        """Generate the table headline."""
        head_line = "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\t"
        if self.rpkm and self.normalization_factors:
            head_line += self._rpkm_column_headers()
        elif self.normalization_factors:
            head_line += self._normalized_couting_column_headers()
        else:
            head_line += self._raw_counting_column_headers()
        return(head_line + "\n")

    def _rpkm_column_headers(self):
        """Generate column header for RPKM values."""
        return("\t".join(
                ["%s RPKM normalized (%s reads)" % (
                        file_name.split("/")[-1], normalization_factor)
                 for file_name, normalization_factor 
                 in zip(self.annotation_mapping_files, 
                        self.normalization_factors)]))

    def _normalized_couting_column_headers(self):
        """Generate column headers for normalized countings."""
        return("\t".join(
                ["%s normalized by %s" % (
                        file_name.split("/")[-1], normalization_factor) 
                 for file_name, normalization_factor 
                 in zip(self.annotation_mapping_files, 
                        self.normalization_factors)]))

    def _raw_counting_column_headers(self):
        """Generate column headers for raw countings."""
        return("\t".join([file_name.split("/")[-1] for 
                   file_name in sorted(self.annotation_mapping_files)]))

    def _calc_rpkm(self, counting, normalization_factor, gene_length):
        """Calculate an RPKM value."""
        if normalization_factor == 0:
            sys.stderr.write(
                "The normalization factor is zero. This should not occure.")
            return(0)
        if gene_length == 0:
            sys.stderr.write(
                "The gene length is zero. This should not occure.")
            return(0)
        return((float(10**9) * float(counting)) / 
               (float(normalization_factor) * float(gene_length)))

    def _read_annotation_mapping_file(self, mapping_file):
        """Read an read mapping annotation overlap file."""
        for line in open(mapping_file):
            if line[0] in ["\n", '#']:
                continue
            entry = self._parse_annotation_mapping_line(line)
            if self._overlap_sufficient(entry) and self._correct_strand(entry):
                self._add_counting(mapping_file, entry)

    def _add_counting(self, mapping_file, entry):
        """Add a overlap counting.

        Mappings of reads that are mapped to more than one place are
        normalized by the total number of mappings of that read.
        """
        add_value = 1.0 # i.e. one read
        if self.count_nucleotides: 
            # The number of overlapping nucleotides
            add_value = self._calc_overlap(entry)
        annotation_entry_key = self._annotation_entry_key(entry)
        self.mapping_files_and_annotation_counting[mapping_file].setdefault(
            annotation_entry_key, 0)
        # The added value is divided by the number of mapping of a read
        # to prevent a too strong influence of single reads.
        # It is also divided by the number of overlap one single mapping
        # has.
        add_value = (
            float(add_value) / 
            float(entry['query_no_of_mappings']) /
            float(entry['no_of_overlaps_of_the_mapping']))
        self.mapping_files_and_annotation_counting[mapping_file][
            annotation_entry_key] += add_value

    def _annotation_entry_key(self, entry):
        """Generate key to descriminate entries."""
        return("%s:%s:%s:%s" % (
            entry["annotation_start"], entry["annotation_end"], 
            entry["annotation_pid"], entry["annotation_gene_name"]))

    def _parse_annotation_mapping_line(self, line):
        """Parse a read mapping annotation overlap line."""
        split_line = line[:-1].split("\t")
        start, end = sorted([int(split_line[1]), int(split_line[2])])
        number_of_mappings = split_line[4].split(":")[-1]
        return({'query_id' : split_line[0],
                'query_start' : start,
                'query_end' : end,
                'query_strand' : split_line[3],
                'query_no_of_mappings' : number_of_mappings,
                'annotation_pid' : split_line[5],
                'annotation_start' : int(split_line[6]),
                'annotation_end' : int(split_line[7]),
                'annotation_strand' : split_line[8],
                'annotation_description' : split_line[9],
                'annotation_gene_name' : split_line[10],
                'annotation_synonym' : split_line[11],
                'annotation_code' : split_line[12],
                'annotation_cog' : split_line[13],
                'no_of_overlaps_of_the_mapping' : split_line[14]
                })

    def _parse_annotation_line(self, line):
        """
        Extracts the information of a annotation line and turn it into
        an entry dictionary.
        """
        split_line = line[:-1].split()
        start, end = sorted([int(pos) for pos in split_line[0].split("..")])
        entry = {
            'annotation_start' : int(start),
            'annotation_end' : int(end),
            'annotation_strand' : split_line[1],
            'annotation_pid' : split_line[3],
            'annotation_gene_name' : split_line[4],
            'annotation_synonym' : split_line[5],
            'annotation_code' : split_line[6],
            'annotation_cog' : split_line[7],
            'annotation_description' : " ". join(split_line[8:])
            }
        return(entry)

    def _overlap_sufficient(self, entry):
        """Test if an overlap is sufficient."""
        overlap = self._calc_overlap(entry)
        if overlap < self.min_overlap:
            return(False)
        if self.min_read_percentage:
            if (self._read_overlap_percentage(entry) 
                < self.min_read_percentage):
                return(False)
        elif overlap < 0:
            sys.stdout.write('Error! No overlap for entry %s.\n' % entry)
            sys.exit(2)
        return(True)

    def _read_overlap_percentage(self, entry):
        """Caculate length percentage of the read mapping overlapping 
        with the annoation.
        """
        return(float(self._calc_overlap(entry)) 
               / float(self._read_mapping_length(entry))
               * 100.0)
    
    def _read_mapping_length(self, entry):
        """Calc the read mapping length."""
        return(int(entry['query_end']) - int(entry['query_start']))

    def _correct_strand(self, entry):
        """Test if the strands have the requested features."""
        if self.strand_orientation == "b":
            return(True)
        elif (self.strand_orientation == "s" and 
            entry["query_strand"] == entry["annotation_strand"]):
            return(True)
        elif (self.strand_orientation == "a" and 
            entry["query_strand"] != entry["annotation_strand"]):
            return(True)
        return(False)
    
    def _calc_overlap(self, entry):
        """Calculate the overlap of a read mapping and an annotation."""
        query_start = entry['query_start']
        query_end = entry['query_end']
        # Start and end setting by min/max is needed as in some file
        # the coordinates of elements on the minus strand are swapped.
        annotation_start = min(entry['annotation_start'],entry['annotation_end'])
        annotation_end = max(entry['annotation_start'],entry['annotation_end'])
        return(min([query_end, annotation_end]) - 
               max([query_start, annotation_start]))

    def _gene_length(self, entry):
        """Calculate the gene length."""
        return(entry['annotation_end'] - entry['annotation_start'])

if __name__ == '__main__': 
    main()
