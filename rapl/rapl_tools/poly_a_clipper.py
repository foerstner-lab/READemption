#!/usr/bin/env python

# Copyright (c) 2011, Konrad Foerstner <konrad@foerstner.org>
#
# Permission to use, copy, modify, and/or distribute this software for
# any purpose with or without fee is hereby granted, provided that the
# above copyright notice and this permission notice appear in all
# copies.

# THE SOFTWARE IS PROVIDED 'AS IS' AND THE AUTHOR DISCLAIMS ALL
# WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
# AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL
# DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
# PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
# TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.

"""Removes the poly A tail of fasta sequences.

usage: poly_a_clipper.py [-h] [-o OUTPUT_FILE_PREFIX] input_fasta_file

positional arguments:
  input_fasta_file

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_FILE_PREFIX, --output_file_prefix OUTPUT_FILE_PREFIX
                        set the name of the output file. If not set theoutput
                        file name is the input file name plus the appendix
                        ".clipped.fa".

"""
__description__ = "Removes the poly A tail of fasta sequences."
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2011 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = "0.1"

import re
import sys
from argparse import ArgumentParser
sys.path.append('../../libs/')
from fasta import FastaParser

def main():
    """Run, forest, run!"""
    parser = create_arg_parser()
    args = parser.parse_args()
    poly_a_clipper = PolyAClipper(
        args.input_fasta_file, args.output_file_prefix)
    poly_a_clipper.read_input_and_clip()

def create_arg_parser():
    parser = ArgumentParser(description=__description__)
    parser.add_argument("input_fasta_file")
    parser.add_argument(
        "-o", "--output_file_prefix", default=None,
        help="set the name of the output file. If not set the" +
        "output file name is the input file name plus the " +
        "appendix \".clipped.fa\".")
    return(parser)

class PolyAClipper(object):
    """Removes the poly A tail of fasta sequences.
    
    [M,D,I]
    ( (AAAAA 1...15 AAA[1,0,0] (G | 0...0) CGGGGCGATGTCTCG[1,1,2] | AAAAA AAAAAAAAA[1,0,0] ) | ( AAAAA AAAAAAAAA[0,0,1] | AAAAA AAAAAAAAA[1,1,0]))

    AAAAA AAAAAAAAA[1,0,0]  => AAAAAAAAA then Hamming distance = 1
    AAAAA AAAAAAAAA[0,0,1]  => AAAAAAAAA + A then Hamming distance = 1?
    AAAAA AAAAAAAAA[1,1,0]  => AAAAAAAAA - A then Hamming distance = 1?

    """

    def __init__(self, input_file, output_file):
        self.input_file = input_file
        if output_file:
            self.output_file = output_file + ".clipped.fa"
        else:
            self.output_file = self.input_file + ".clipped.fa"
     
    def read_input_and_clip(self):
        """Read entry by entry and returns the clipped version.

        Write output on the fly.
        """
        fasta_parser = FastaParser()
        output_fh = open(self.output_file, "w")
        for header, seq in fasta_parser.parse_fasta_file(self.input_file):
            output_fh.write(">%s\n%s\n" % (header, self._test_and_clip(seq)))
        output_fh.close()
        sys.stdout.write("Wrote output to: %s\n" % self.output_file)

    def _test_and_clip(self, sequence):
        """Return poly A free sequences.

        Clip the sequence if poly-A tail exists - return the
        sequence unproccessed if this is not the case.
        """
        sequence = sequence.upper()
        length = 11
        if "AAAA" in sequence:
            for subseq, start_pos in self._aaaa_starting_substrings(
                sequence, length):
                # Tolerate one mismatch
                if subseq.count("A") < length - 1:
                    continue
                else:
                    # Use sequence only to the start of the poly-A
                    # tail
                    sequence = sequence[:start_pos]
                    # If the poly A is at the start of a sequence set
                    # the sequence to a single A
                    if start_pos == 0:
                        sequence = "A"
                    break
        sequence = self._remove_3_prime_a(sequence)
        return(sequence)

    def _aaaa_starting_substrings(self, sequence, length):
        start_pos = 0
        while start_pos != -1:
            start_pos = sequence.find("AAAA", start_pos)
            if start_pos != -1:
                if start_pos + length > len(sequence):
                    start_pos = -1
                else:
                    cur_start_pos = start_pos
                    start_pos = start_pos + 1
                    yield([sequence[cur_start_pos:cur_start_pos+length], 
                           cur_start_pos])

    def _remove_3_prime_a(self, sequence):
        if sequence == '':
            return(sequence)
        elif sequence[-1] == "A":
            sequence = self._remove_3_prime_a(sequence[:-1])
        return(sequence)
                    
if __name__ == '__main__': 
    main()
