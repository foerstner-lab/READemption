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
Usage: filter_fasta_entries_by_size.py [options] <FASTA FILE> <MIN SEQ LENGTH>

Filters fasta sequences by size.

Options:
  -h, --help            show this help message and exit
  -o OUTPUT_FILE_PREFIX, --output_file_prefix=OUTPUT_FILE_PREFIX
         
"""
__description__ = "Filters fasta sequences by size."
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2011 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = "0.1"

import sys
from optparse import OptionParser
sys.path.append('../../libs/')
from fasta import FastaParser

def main():
    """Run, forest, run!"""
    python_version = sys.version.split()[0]
    if not python_version[0] == "3":
        sys.stdout.write('Error! Please us Python 3.x!\n')
        sys.exit(2)
    parser = OptionParser(
        "Usage: %prog [options] <FASTA FILE> <MIN SEQ LENGTH>\n\n" + 
        __description__)
    parser.add_option("-o", "--output_file_prefix", default=None, 
                      action="store", dest="output_file_prefix", help="")
    (options, args) = parser.parse_args()
    if not len(args) == 2:
        parser.print_help()
        sys.exit(2)
    faster_size_filter = FasterSizeFilter(args, options)
    faster_size_filter.read_fasta_file_and_print()

class FasterSizeFilter(object):
    """Converts segemehl output to FASTA format.

    Two files are generated. One that contain entries that are longer
    than or equal the given cut-off, another one with entries which do
    not fullfil the length criterion.
    
    """

    def __init__(self, args, options):
        """Create instance.

        Arguments:
        - args:
          - args[0]: the input fasta file name
          - args[1]: mininmal sequence length
        - options:
          - options.output_file_prefix: Set a output file prefix. If not 
                                        given, the input file path is used to
                                        generate the output file path.

        """
        self.fasta_file = args[0]
        self.min_nucleotides = int(args[1])
        if options.output_file_prefix:
            output_file_start = options.output_file_prefix
        else:
            output_file_start = args[0]
        self.output_file_sufficent = "%s.size_filtered_gtoe_%sbp.fa" % (
            output_file_start, self.min_nucleotides)
        self.output_file_insufficent = "%s.size_filtered_lt_%sbp.fa" % (
            output_file_start, self.min_nucleotides)

    def read_fasta_file_and_print(self):
        """Read the input file and do the filtering on the fly. """
        fasta_parser = FastaParser()
        output_fh_sufficient = open(self.output_file_sufficent, "w")
        output_fh_insufficient = open(self.output_file_insufficent, "w")
        counter_sufficient = 0
        counter_insufficient = 0
        for header, seq in fasta_parser.parse_fasta_file(self.fasta_file):
            if len(seq) >= self.min_nucleotides:
                output_fh_sufficient.write(">%s\n%s\n" % (header, seq))
                counter_sufficient += 1
            else:
                output_fh_insufficient.write(">%s\n%s\n" % (header, seq))
                counter_insufficient += 1
        output_fh_sufficient.close()
        output_fh_insufficient.close()
        sys.stdout.write(
            "Wrote %s entries to file \"%s\" (>= %s bp).\n" 
            "Wrote %s entries to file \"%s\" (< %s bp).\n" % (
                counter_sufficient, self.output_file_sufficent, 
                self.min_nucleotides, counter_insufficient, 
                self.output_file_insufficent, self.min_nucleotides))
                
if __name__ == '__main__': 
    main()
