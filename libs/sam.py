"""
Copyright (c) 2011, Konrad Foerstner <konrad@foerstner.org>

Permission to use, copy, modify, and/or distribute this software for
any purpose with or without fee is hereby granted, provided that the
above copyright notice and this permission notice appear in all
copies.

THE SOFTWARE IS PROVIDED 'AS IS' AND THE AUTHOR DISCLAIMS ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL
DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.
         
"""
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2011 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"

import re

class SamParser(object):
    """ A class for parsing and handling SAM files produced by segemehl. 

    Really limited to the output of segemehl only!
    """

    def entries(self, sam_file):
        for line in self._entry_lines(sam_file):
            yield(self._line_to_entry(line))

    def _entry_lines(self, sam_file):
        for line in open(sam_file):
            if line[0] in ["\n", "@"]:
                continue
            yield(line)

    def _line_to_entry(self, line):
        line = line.strip()
        split_line = line.split("\t")
        entry = {
            'query' : split_line[0],
            'flag' : int(split_line[1]),
            'reference' : split_line[2],
            'start' : int(split_line[3]),
            'mapping_quality' : int(split_line[4]),
            'cigar' : split_line[5],
            'mate' : split_line[6],
            'mate_start' : int(split_line[7]),
            'template_length' : int(split_line[8]),
            'sequence' : split_line[9],
            'phred_quality' : split_line[10],
            'distance' : split_line[11],
            'mismatches' : split_line[12],
            'number_of_hits' : split_line[13]
            }
        return(entry)

    def number_of_hits_as_int(self, entry):
        return(int(entry['number_of_hits'].split(":")[-1]))

    def entry_start_end_strand(self, entry):
        try:
            self.m_pattern
            self.d_pattern
        except AttributeError:
            self.m_pattern = re.compile("(\d+)M")
            self.d_pattern = re.compile("(\d+)D")
        matches_mismatches = [int(x) for x in re.findall(
                self.m_pattern, entry["cigar"])]
        deletions = [int(x) for x in re.findall(
                self.d_pattern, entry["cigar"])]
        end = entry["start"] + sum(matches_mismatches) + sum(deletions) - 1
        return(entry["start"], end, self._flag_to_strand(entry["flag"]))
    
    def _flag_to_strand(self, flag):
        if flag == 0:
            return("+")
        elif flag == 16:
            return("-")
        else:
            return("ERROR")

    def header_lines(self, sam_file):
        head = ""
        for line in open(sam_file):
            if line.startswith("\n"):
                continue
            elif line.startswith("@"):
                head += line
            else:
                exit
        return(head)
    
    def reference_list(self, sam_file):
        """Return a list of all reference sequence in a SAM file."""
        references = {}
        for entry in self.entries(sam_file):
            references[entry['reference']] = 1
        return(sorted(references.keys()))

class SamBuilder(object):

    def entry_to_line(self, entry):
        return("%s\n" % (
                "\t".join(
                    [entry["query"],
                     str(entry["flag"]),
                     entry["reference"],
                     str(entry["start"]),
                     str(entry["mapping_quality"]),
                     entry["cigar"],
                     str(entry["mate"]),
                     str(entry["mate_start"]),
                     str(entry["template_length"]),
                     entry["sequence"],
                     entry["phred_quality"],
                     entry["distance"],
                     entry["mismatches"],
                     entry["number_of_hits"]])))
               
