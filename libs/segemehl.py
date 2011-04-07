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

class SegemehlParser:
    """ A class for parsing and handling segemehl files. """

    def entries(self, segemehl_file):
        for line in self._entry_lines(segemehl_file):
            yield(self._line_to_entry(line))

    def _entry_lines(self, segemehl_file):
        for line in open(segemehl_file):
            if line[0] in ["\n", "#"]:
                continue
            yield(line)

    def _line_to_entry(self, line):
        line.strip()
        split_line = line.split("\t")
        entry = {
            'id' : split_line[0],
            'semiglobal_alignment_edist' : int(split_line[1]),
            'seed_score' : int(split_line[2]),
            'seed_evalue' : float(split_line[3]),
            'seed_start' : int(split_line[4]),
            'seed_end' : int(split_line[5]),
            'matches' : int(split_line[6]),
            'mismatches' : int(split_line[7]),
            'insertions' : int(split_line[8]),
            'deletions' : int(split_line[9]),
            'strand' : split_line[10],
            'hit_start' : int(split_line[11]),
            'hit_end' : int(split_line[12]),
            'target_description' : split_line[13],
            'matching_code' : split_line[14],
            'sequence' : split_line[15].strip()
            }
        return(entry)

    def header(self):
        pass

class SegemehlBuilder:
    ""
    
    def entry_to_line(self, entry):
        return("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t"
               "%s\t\n" % (
                entry['id'],
                entry['semiglobal_alignment_edist'],
                entry['seed_score'],
                entry['seed_evalue'],
                entry['seed_start'],
                entry['seed_end'],
                entry['matches'],
                entry['mismatches'],
                entry['insertions'],
                entry['deletions'],
                entry['strand'],
                entry['hit_start'],
                entry['hit_end'],
                entry['target_description'],
                entry['matching_code'],
                entry['sequence']))



