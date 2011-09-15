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
# DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
# PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
# TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.

"""
Filters SAM entries by the percentage of a nucleotide.

Usage: filter_segemehl_by_nucleotide_percentage.py [options] <SEGEMEHL FILE> <NUCLEOTIDE> <PERCENTAGE>

Options:
  -h, --help            show this help message and exit
  -o OUTPUT_FILE_PREFIX, --output_file_prefix=OUTPUT_FILE_PREFIX

"""
__description__ = "Filters Segemehl entries by the percentage of a nucleotide."
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2011 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = "0.1"

import sys
from optparse import OptionParser
sys.path.append('../libs')
from sam import SamParser
from sam import SamBuilder

def main():
    """Run, forest, run!"""
    python_version = sys.version.split()[0] 
    if not python_version[0] == "3":
        sys.stdout.write('Error! Please us Python 3.x!\n')
        sys.exit(2)
    parser = OptionParser(
        "Usage: %prog [options] <SEGEMEHL FILE> <NUCLEOTIDE> <PERCENTAGE>\n\n" 
        + __description__ )
    parser.add_option("-o", "--output_file_prefix", default=None, action="store", 
                      dest="output_file_prefix", help="output file prefix. If "
                      "not set the input file path will be used instead.")
    (options, args) = parser.parse_args()
    if not len(args) == 3:
        parser.print_help()
        sys.exit(2)
    sam_nucl_perc_filter = SamNuclPercFilter(args, options)
    sam_nucl_perc_filter.read_mapping_file_and_print()

class SamNuclPercFilter(object):
    """Filters SAM entries by the percentage of a nucleotide.

    It uses the complementary base when analysing read from the minus
    strand.
    """
 
    def __init__(self, args, options):
        """Create instance.

        Arguments:
        - args:
          - args[0]: SAM mapping file
          - args[1]: The nucleotide that is used for the percentage
                     caculation (A, C, G or T)
          - args[2]: The            
        - options: 
          - options.output_file_prefix: The prefix for the output files.
        """
        
        self.mapping_file = args[0]
        self.nucleotide = args[1].upper()
        if self.nucleotide not in ['A', 'C', 'G', 'T']:
            sys.stderr.write(
                "Error: %s is not a valid nucleotide.\n" % self.nucleotide)
            sys.exit(2)
        self.max_percentage = float(args[2])
        if (self.max_percentage < 0 or self.max_percentage > 100):
            sys.stderr.write(
                "Error: The maximum percentage must be between 0 and 100.\n")
            sys.exit(2)
        if options.output_file_prefix:
            output_file_start = options.output_file_prefix
        else:
            output_file_start = args[0]
        # ltoe = less than or equal
        self.output_file_sufficent = "%s.filtered_ltoe_%s%%_%s.txt" % (
            output_file_start, self.max_percentage, self.nucleotide)
        # gt = greater than
        self.output_file_insufficent = "%s.filtered_gt_%s%%_%s.txt" % (
            output_file_start, self.max_percentage, self.nucleotide)
        self.sam_parser = SamParser()
        self.sam_builder = SamBuilder()
        self.counter_sufficient = 0
        self.counter_insufficient = 0

    def read_mapping_file_and_print(self):
        """Read the mapping file and write the output on the fly.

        Two output files are generated. One for the entries that have
        the maximum percentage or below, another on for the one that
        have percentage values above the given cut-off.

        """
        self._write_parameter_settings()
        self._setup_output_files()
        for entry in self.sam_parser.entries(self.mapping_file):
            self._process_sam_entry(entry)
        self.output_fh_sufficient.close()
        self.output_fh_insufficient.close()
        self._write_processing_stats()
        
    def _write_parameter_settings(self):
        sys.stdout.write("Reading file: %s\n" % self.mapping_file)
        sys.stdout.write("Searching the nucleotide: %s\n" % self.nucleotide)
        sys.stdout.write("Max percentage: %s\n" % self.max_percentage)        

    def _setup_output_files(self):
        self.output_fh_sufficient = open(self.output_file_sufficent, "w")
        self.output_fh_insufficient = open(self.output_file_insufficent, "w")
        sam_header = self.sam_parser.header_lines(self.mapping_file)
        self.output_fh_sufficient.write(sam_header)
        self.output_fh_insufficient.write(sam_header)

    def _process_sam_entry(self, entry):
        line = self.sam_builder.entry_to_line(entry)
        if self._nucl_percentage_too_high(
            entry["sequence"], self.sam_parser._flag_to_strand(entry["flag"])):
            self.output_fh_insufficient.write(line)
            self.counter_insufficient += 1
        else:
            self.output_fh_sufficient.write(line)
            self.counter_sufficient += 1

    def _write_processing_stats(self):
        sys.stdout.write(
            "Wrote %s entries to file \"%s\" (<= %s %% ).\n" 
            "Wrote %s entries to file \"%s\" (> %s %%).\n" % (
                self.counter_sufficient, self.output_file_sufficent, 
                self.max_percentage, self.counter_insufficient, 
                self.output_file_insufficent, self.max_percentage))        

    def _nucl_percentage_too_high(self, sequence, strand):
        """Test if the percentage exceed the cut-off value."""
        sequence = sequence.upper()
        nucleotide = self._effective_nucleotide(strand)
        return(self._perc_in_sequence(nucleotide, sequence) > self.max_percentage)

    def _effective_nucleotide(self, strand):
        """Set the character to search.

        Takes the complementary nucleotide if the minus strand is
        scanned.
        """
        if strand == "+":
            return(self.nucleotide)
        elif strand == "-":
            return({"A" : "T", "C" : "G", "G" : "C", "T" : "A"}[
                    self.nucleotide])

    def _perc_in_sequence(self, nucleotide, sequence):
        return(float(sequence.count(nucleotide) / float(len(sequence)) * 100))
    
if __name__ == '__main__': 
    main()
