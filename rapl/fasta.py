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

class FastaParser:
    """ A class for parsing and handling fasta files. """

    def parse_fasta_file(self, fasta_file_name):
        """
        Go line by line though the file and return fasta
        entries. Generator based implementation, Memory saving
        returns header and sequence for each iterator.
        """
        current_header = ''
        current_sequence = ''

        # Switch telling if this the first entry
        virgin = True 
        for line in open(fasta_file_name):
            # For usual headers
            if line[0] == '>' and virgin == False:
                yield current_header, current_sequence
                current_header = line[1:-1]
                current_sequence = ''
            # For the very first header of the file
            elif line[0] == '>' and virgin == True:
                virgin = False
                current_header = line[1:-1]
            # For sequence lines
            else:
                current_sequence += line[:-1]
        # For the last entry
        yield current_header, current_sequence

    def shape_sequence(self, sequence, length):
        """ Breaks the sequence into lines of a given length."""
        counter = 0
        shaped_sequence = ''
        for character in sequence:
            counter += 1
            shaped_sequence += character
            if counter == length:
                shaped_sequence += '\n'
                counter = 0
        if not len(shaped_sequence) == 0:
            if shaped_sequence[-1] == '\n':
                shaped_sequence = shaped_sequence[:-1]                
        else:
            print('Warning: Sequence lenght = 0!')
        return(shaped_sequence)
