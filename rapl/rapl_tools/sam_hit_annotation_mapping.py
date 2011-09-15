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

"""Find overlaps of read mappings and annotation table entries
(e.g. an ptt. or .rna file).

usage: sam_hit_annotation_mapping.py [-h] [-m MIN_OVERLAP] [-o OUTPUT_FILE]
                                     [-v]
                                     SAM_FILE ANNOTATION_FILE REFERENCE_HEADER

Find overlaps of read mappings and annotation table entries (e.g. an ptt. or
.rna file).

positional arguments:
  SAM_FILE              a read mapping file (e.g. from segemehl)
  ANNOTATION_FILE       a (NCBI) annotation file
  REFERENCE_HEADER      the header of the reference genome file that should be
                        considered. Only mapping to that genome file will be
                        placed to annotations.

optional arguments:
  -h, --help            show this help message and exit
  -m MIN_OVERLAP, --min_overlap MIN_OVERLAP
                        minimal overlap needed to count a match (in bp).
                        Default is 1.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        set the name of the output file. If not set the output
                        file name is the annotaion file name plus the appendix
                        "_annotation_mapping".
  -v, --version         show program's version number and exit
"""

__description__ = ("Find overlaps of read mappings and annotation " +
                   "table entries (e.g. an ptt. or .rna file).")
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2011 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = "0.2"

import bisect
import sys
from argparse import ArgumentParser
sys.path.append('../libs')
from sam import SamParser

def main():
    """Run, forest, run!"""
    python_version = sys.version.split()[0]
    if not python_version[0] == "3":
        sys.stdout.write('Error! Please us Python 3.x!\n')
        sys.exit(2)
    parser = ArgumentParser(description=__description__)
    parser.add_argument("sam_file", metavar='SAM_FILE', 
                        type=str, 
                        help="a read mapping file (e.g. from segemehl)")
    parser.add_argument("annotation_file", metavar='ANNOTATION_FILE', 
                        type=str, help="a (NCBI) annotation file")
    parser.add_argument("reference_header", metavar='REFERENCE_HEADER', 
                        type=str, help="the header of the reference genome "
                        "file that should be considered. Only mapping to that "
                        "genome file will be placed to annotations.")
    parser.add_argument("-m", "--min_overlap", default=1, action="store",
                        type=int, dest="min_overlap",  help="minimal overlap "+
                        "needed to count a match (in bp). Default is 1.")
    parser.add_argument("-o", "--output_file", action="store", 
                        default=None, dest="output_file", type=str,
                        help="set the name of the output file. If not set the " +
                        "output file name is the annotaion file name plus the "+
                        "appendix \"_annotation_mapping\".")
    parser.add_argument("-v", "--version", action='version', version=__version__)
    args = parser.parse_args()
    sam_hit_annotation_mapping = SamHitAnnotationMapping(
        args.sam_file, args.annotation_file, args.reference_header, 
        args.output_file,  args.min_overlap)
    sam_hit_annotation_mapping.map()

class SamHitAnnotationMapping(object):
    """Maps sam hits to an annotation table entries (in e.g. an
    ptt. or .rnt file).

    Overlaps of sam hits and annotation entries are searched.
    """

    def __init__(
        self, sam_file, annotation_file, reference_header, output_file, 
        min_overlap):
        """Create instance.

        Arguments:
        - sam_file
        - annotation_file

        """
        self.sam_file = sam_file
        self.annotation_file = annotation_file
        self.reference_header = reference_header
        self.min_overlap = int(min_overlap)
        self.output_file = output_file
        self.annotations = []
        self.sorted_annotations = []
        self.sorted_annotation_starts = []
        
        # Two parameter to guarantee the catching of all annotation
        # hits.
        self.steps_back = 5
        self.max_number_of_non_overlappings = self.steps_back + 5

    def map(self):
        """Do the actual mapping."""
        self._get_annotation_entries()
        self._sort_annotation_entries_by_start()
        self._generate_output_file()
        sam_parser = SamParser()
        for raw_sam_entry in sam_parser.entries(self.sam_file):
            # Skip entries that represent mappings to other reference
            # genome files
            if not raw_sam_entry["reference"] == self.reference_header:
                continue
            start, end, strand = sam_parser.entry_start_end_strand(raw_sam_entry)
            number_of_hits = raw_sam_entry['number_of_hits']
            sam_entry = {
                "start" : start, "end" : end, "strand" : strand,  
                "id" : raw_sam_entry["query"],      'number_of_hits' : number_of_hits}
            for hit_annotation_entry in self._search_in_list(sam_entry):
                self._write_to_output_file(sam_entry, hit_annotation_entry)
    def _get_annotation_entries(self):
        """Collect the annotation entries"""

        for entry in self._parse_annotation_file(self.annotation_file):
            self.annotations.append(entry)

    def _sort_annotation_entries_by_start(self):
        self.sorted_annotations = sorted(
            self.annotations, key=lambda k: [k['start'], k["end"]])
        self.sorted_annotation_starts = [
            annotation["start"] for annotation in self.sorted_annotations]

    def _parse_annotation_file(self, input_file):
        """ Reads annotation table file and return processed entries."""
        for line in open(input_file):
            split_line = line.split("\t")
            if '..' not in split_line[0] or len(split_line) != 9:
                continue
            yield(self._table_line_to_entry(line))

    def _table_line_to_entry(self, line):
        """Extract the information of a annotation line and turn it
        into an entry dictionary.
        
        """
        split_line = line[:-1].split()
        start, end = sorted(int(pos) for pos in split_line[0].split(".."))
        entry = {'start' : start, 'end' : end, 'strand' : split_line[1],
            'pid' : split_line[3], 'gene_name' : split_line[4],
            'synonym' : split_line[5], 'code' : split_line[6],
            'cog' : split_line[7], 'description' : " ". join(split_line[8:])}
        return(entry)

    def _search_in_list(self, sam_hit):
        """Searches for overlaps of all annotation entries in a given
        list with a given sam entry.
        
        """
        index = self._annotation_index(sam_hit)
        number_of_non_overlappings = 0
        for annotation_entry in self.sorted_annotations[index:]:
            if self._has_sufficiant_overlap(sam_hit, annotation_entry):
                number_of_non_overlappings = 0
                yield(annotation_entry)
            else:
                # Go a certain amount of annotations without overlap
                # ahead to be sure not to miss one.
                number_of_non_overlappings += 1
                if number_of_non_overlappings > self.max_number_of_non_overlappings:
                    break

    def _annotation_index(self, sam_hit):
        # The index of the annotation with the start position that is
        # closest (left side. i.e. larger) to the start position of the
        # read which should be checked.
        index = bisect.bisect_left(
            self.sorted_annotation_starts, sam_hit["start"])
        if index == len(self.sorted_annotation_starts):
            index -= 1
        # For a set of annotations that start with the same position
        while (index > 0 and self.sorted_annotation_starts[index] == 
               self.sorted_annotation_starts[index-1]):
            index -= 1
        # The steps back are done to cover the constellation in which
        # the read start position value is larger than the starting
        # position of the overlapping annoation.
        # 
        #  -------  annotation
        #     -------      read
        if index - self.steps_back > 0:
            index = index - self.steps_back
        else:
            index = 0
        return(index)

    def _write_to_output_file(self, sam_hit, annotation_entry):
        """Wrtie the results to the output file."""
        self.output_fh.write(
            "\t".join([str(sam_hit[key]) for key in [
                        "id", "start", "end", "strand", "number_of_hits"]])
            + "\t" 
            + "\t".join([str(annotation_entry[key]) for key in [
                        "pid", "start", "end", "strand", "description",
                        "gene_name", "synonym", "code", "cog"]])
            + "\n")

    def _generate_output_file(self):
        """Set up an output file and add a header line."""
        if self.output_file:
            self.output_fh = open(self.output_file, "w")
        else:
            self.output_fh = open("%s.annotation_mapping" % 
                                  self.annotation_file, "w")
        self.output_fh.write(
            "\t".join([
                "#Sam hit id", "Sam hit start", "Sam hit end", "Sam hit strand",
                "Sam number of hits of read", "Annotation PID", 
                "Annotation start", "Annotation end", "Annotation strand",
                "Annotation description", "Annotation Gene name",
                "Annotation Synonym", "Annotation Code", "Annotation COG"]) 
            + "\n")

    def _add_annotation_counting(self, description, is_protein):
        """Increaeses the counter of the group of a given description
        by one.

        """
        if is_protein:
            description = 'protein'
        else:
            if description.startswith('Anticodon:') or "tRNA" in description:
                description = 'tRNA'
        self.annotation_counting.setdefault(description, 0)
        self.annotation_counting[description] += 1

    def _has_sufficiant_overlap(self, hit_entry, annotation_entry):
        """Test if a sam hit and an annotation entry overlap."""
        # Start and end setting by min/max is needed as in some file
        # the coordinates of elements on the minus strand are swapped.
        annotation_start = min(annotation_entry['start'], 
                               annotation_entry['end'])
        annotation_end = max(annotation_entry['start'], 
                             annotation_entry['end'])
        ## Overlap must be greater than 0
        #return(0 <= (min([hit_entry['end'], annotation_end]) - 
        #       max([hit_entry['start'], annotation_start])))        
        # Overlap must be equal or greater than the given cut-off value
        return(self.min_overlap <= (
                min([hit_entry['end'], annotation_end]) - 
                max([hit_entry['start'], annotation_start])))        

if __name__ == '__main__': 
    main()
