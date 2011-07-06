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

"""Converts segemehl SAM output to GR format suited for the Integrated
Genome Browser.

usage: sam2gr.py [-h] [-t MAPPING_TARGET] [-r] [-n] [-f FLOAT] [-m FLOAT]
                 [-o PATH] [-s]
                 SAM_FILE

positional arguments:
  SAM_FILE              A SAM file.

optional arguments:
  -h, --help            show this help message and exit
  -t MAPPING_TARGET, --target MAPPING_TARGET
                        If the SAM file contains mapping to more than one
                        genetic element (chromosomes, plasmids etc.) one
                        target/reference element has to be seleceted by givin
                        the header string of it.
  -r, --normalize_by_reads
                        normalize each intensity by number of reads (not
                        mappings!) and multiply by 1000. Thus the y-values are
                        per mill. The multiplier can be specified using -m.
  -n, --normalize_by_nucleotides
                        normalize each intensity by number of mapped
                        nucleotides and multiply by 1000. Thus the y-values
                        are per mill. The multiplier can be specified using
                        -m.
  -f FLOAT, --normalization_factor FLOAT
                        normalized with given factor.
  -m FLOAT, --multiplier FLOAT
                        set the multiplier for the normalization. If not
                        specified the normalize_by_readsd value is multiplied
                        by 1000. Only used if -n or -r is set.
  -o PATH, --output_prefix PATH
                        output path/name prefix. If not specified the input
                        file/path will be used and extended by an appendix.
  -s, --starts_only     Count only the first base of a read mapping. This can
                        be as simple approach to detect transcription start
                        sites (TSS).
"""

__description__ = ("Converts segemehl output to GR format suited for the " +
                   "Integrated Genome Browser.")
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2011 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = "0.1"

import sys
from argparse import ArgumentParser
from sam import SamParser

def main():
    """Run, forest, run!"""
    python_version = sys.version.split()[0]
    if not python_version[0] == "3":
        sys.stdout.write('Error! Please us Python 3.x!\n')
        sys.exit(2)
    arg_parser = set_args()
    args = arg_parser.parse_args()
    sam2gr = Sam2Gr(
        args.sam_file, args.mapping_target, args.normalize_by_reads,
        args.normalize_by_nucleotides, args.normalization_factor,
        args.normalization_multiplier, args.output_prefix, 
        args.starts_only)
    sam2gr.check_parameters()
    sam2gr.collect_intensities()
    sam2gr.print_output()

def set_args():
    parser = ArgumentParser(description=__description__)
    parser.add_argument("sam_file", metavar='SAM_FILE', help="A SAM file.")
    parser.add_argument("-t", "--target", default=None, dest="mapping_target", 
                        action="store", help="If the SAM file contains mapping " + 
                        "to more than one genetic element (chromosomes, " + 
                        "plasmids etc.) one target/reference element has to be "
                        "seleceted by givin the header string of it.")
    parser.add_argument("-r", "--normalize_by_reads", default=False, 
                        dest="normalize_by_reads",
                        help="normalize each intensity by number of reads "
                        "(not mappings!) and multiply by 1000. Thus the "
                        "y-values are per mill. The multiplier can be "
                        "specified using -m.", action="store_true")
    parser.add_argument("-n", "--normalize_by_nucleotides", 
                        default=False, dest="normalize_by_nucleotides",
                        help="normalize each intensity by number of mapped "
                        "nucleotides and multiply by 1000. Thus the "
                        "y-values are per mill. The multiplier can be "
                        "specified using -m.", action="store_true")
    parser.add_argument("-f", "--normalization_factor", default=None, 
                        dest="normalization_factor",
                        help="normalized with given factor.", 
                        metavar="FLOAT")
    parser.add_argument("-m", "--multiplier", default=None, 
                        dest="normalization_multiplier",
                        help="set the multiplier for the normalization. "
                        "If not specified the normalize_by_readsd value is " +
                        "multiplied by 1000. Only used if -n or -r is set.", 
                        metavar="FLOAT")
    parser.add_argument("-o", "--output_prefix", default=None, 
                        dest="output_prefix",
                        help="output path/name prefix. If not specified " +
                        "the input file/path will be used and extended by an " +
                        "appendix.", metavar="PATH")
    parser.add_argument("-s", "--starts_only", default=False, 
                        dest="starts_only", action="store_true",
                        help="Count only the first base of a read mapping. " +
                        "This can be as simple approach to detect " +
                        "transcription start sites (TSS).")
    return(parser)

class Sam2Gr(object):
    """Converts sam output to GR format.

    The GR format is suited for the Integrated Genome Browser.
    """

    def __init__(
        self, input_file, mapping_target=None, normalize_by_reads=False, 
        normalize_by_nucleotides=False, normalization_factor=None, 
        normalization_multiplier=None, output_prefix=None, starts_only=False):
        """Create the converter

        Arguments
        
        """
        self.input_file = input_file
        self.mapping_target = mapping_target
        self.normalize_by_reads = normalize_by_reads
        self.normalize_by_nucleotides = normalize_by_nucleotides
        self.normalization_factor = normalization_factor
        self.normalization_multiplier = normalization_multiplier
        self.starts_only = starts_only
        if output_prefix:
            self.output_prefix = output_prefix
        else:
            self.output_prefix = self.input_file
        self.final_normalization_factor = None
        self.intensities_plus = []
        self.intensities_minus = []
        self.reads_and_countings = {}
        self.sam_parser = SamParser()
        self.default_normalization_multiplier = 1000.0

    def check_parameters(self):
        """Check if the given parameters make sense."""
        if (self.normalization_factor != None and 
            float(self.normalization_factor) == float(0)):
            sys.stderr.write(
                "The normalization factor is 0 - this does not make "
                "any sense.\n")
            sys.exit(2)
        if (self.normalization_multiplier != None and 
            float(self.normalization_multiplier) == float(0)):
            sys.stderr.write(
                "The normalization multiplier is 0 - this does not make "
                "any sense.\n")
            sys.exit(2)

    def count_nucleotides(self):
        """Count the number of mapped nucleotides
        
        For each read only one mapping is considered.
        """
        no_of_nucleotides = 0
        prev_entry = ""
        for entry in self.sam_parser.entries(self.input_file):
            if entry['query'] == prev_entry: continue
            no_of_nucleotides += len(entry["sequence"])
            prev_entry = entry['query']
        return(no_of_nucleotides)

    def collect_intensities(self):
        """Read the SAM file and add mapping hits to the coordinates."""
        self._check_if_reference_string_is_needed()
        for entry in self.sam_parser.entries(self.input_file):
            self._process_entry(entry)

    def _process_entry(self, entry):
        # If there are more than one mapping target in the file
        # use only the one with the specified reference.
        if (self.mapping_target and 
            not entry["reference"] == self.mapping_target):
            return()
        start, end, strand = self.sam_parser.entry_start_end_strand(entry)
        number_of_hits = self.sam_parser.number_of_hits_as_int(entry)
        if strand == '+':
            self._add_value_to_coordinates(
                self.intensities_plus, start, end, strand, number_of_hits)
        elif strand == '-':
            self._add_value_to_coordinates(
                self.intensities_minus, start, end, strand, number_of_hits)

    def _check_if_reference_string_is_needed(self):
        if not self.mapping_target:
            references = self.sam_parser.reference_list(self.input_file)
            if len(references) > 1:
                sys.stderr.write(
                    "More than one reference files (plasmid, chromosom) were " 
                    "used to generate this SAM file. Please specify which " 
                    "mappings should be used by using the -t option with the " 
                    "header of the file. The following reference file headers "
                    "were found:\n")
                for reference in references:
                    sys.stderr.write("- %s\n" % reference)
                sys.exit(2)

    def print_output(self):
        """Print two output file - one for each strand."""
        # With normalization
        if self._use_normalization():
            self._set_final_normalization_factor()
            sys.stdout.write("Normalization factor: %s\n" % (
                    self.normalization_factor))
            sys.stdout.write("Normalization multiplier: %s\n" % (
                    self.normalization_multiplier))
        # Without normalization
        else:
            sys.stdout.write(
                "No normalization done. Writing raw intensities.\n")
        suffix_plus = self._generate_suffix("plus")
        self._print_to_file(self.intensities_plus, suffix_plus)
        suffix_minus = self._generate_suffix("minus")
        self._print_to_file(self.intensities_minus, suffix_minus)

    def _add_value_to_coordinates(
        self, intensities, start, end, strand, number_of_hits):
        """Add 1 to each position from start to end position."""
        if len(intensities) < end:
            intensities.extend(
                self._extension_for_intensities_list(intensities, end))
        if self.starts_only:
            if strand == "+":
                intensities[start-1] += 1.0 / float(number_of_hits)
            elif strand == "-":
                intensities[end-1] -= 1.0 / float(number_of_hits)
        else:
            for position in range(start-1, end):
                if strand == "+":
                    intensities[position] += 1.0 / float(number_of_hits)
                elif strand == "-":
                    intensities[position] -= 1.0 / float(number_of_hits)

    def _extension_for_intensities_list(self, intensities, end):
        """Return a list that extends the intensitiy list."""
        return([0] * (end - len(intensities)))

    def _generate_suffix(self, strand):
        """Generate a output file suffix dependend of the normalization."""
        if self._use_normalization() :
            return("_normalized_by_%s_multiplied_by_%s_%s.gr" % (
                    self.normalization_factor_for_suffix, 
                    self.normalization_multiplier_for_suffix, strand))
        else:
            return("_%s.gr" % strand)

    def _print_to_file(self, intensities, file_suffix):
        """Print all intensities to a file with a given suffix. """
        file_name = "%s%s" % (self.output_prefix, file_suffix)
        out_fh = open(file_name, "w")
        position = 1
        # With normalization
        if self._use_normalization():
            for intensity in intensities:
                if intensity != 0.0: # To save space skipt 0 value entries
                    out_fh.write("%s\t%s\n" % (
                            position, self._normalized_intensity(intensity)))
                position += 1
        # Without normalization
        else:
            for intensity in intensities:
                if intensity != 0.0: # To save space skipt 0 value entries
                    out_fh.write("%s\t%s\n" % (position, intensity))
                position += 1
        out_fh.close()
        sys.stdout.write("Wrote to file \"%s\".\n" % file_name)
        
    def _normalized_intensity(self, intensity):
        """Calculate the normalized intensity."""
        return(intensity * self.final_normalization_factor)
    
    def _set_final_normalization_factor(self):
        """Calculate and set the combined normalization factor. """
        # Set the normalization factor
        normalization_factor = None
        if self.normalize_by_nucleotides:
            normalization_factor = self.count_nucleotides()
        elif self.normalize_by_reads:
            normalization_factor = self.sam_parser.number_of_mapped_reads(
                self.input_file)
        elif self.normalization_factor:
            normalization_factor = float(self.normalization_factor)
        # Set the multiplier
        normalization_multiplier = self.default_normalization_multiplier
        if self.normalization_multiplier:
            normalization_multiplier = float(self.normalization_multiplier)
        # Calculate the combined factor
        try:
            self.final_normalization_factor = (
                normalization_multiplier / normalization_factor)
        except ZeroDivisionError:
            sys.stdout.write("Normalization factor 0. Now set to 1.\n")
            self.final_normalization_factor = self.normalization_multiplier
            normalization_factor = 1
        # Needed for suffix creation:
        self.normalization_factor_for_suffix = normalization_factor
        self.normalization_multiplier_for_suffix = normalization_multiplier
        
    def _use_normalization(self):
        return(self.normalize_by_reads or 
               self.normalize_by_nucleotides or
               self.normalization_factor)

if __name__ == '__main__': 
    main()
