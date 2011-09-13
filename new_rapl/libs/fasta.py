class FastaParser(object):
    """ A class for parsing and handling fasta files. """

    def entries(self, fasta_fh):
        """Go line by line though the file and return fasta entries. """
        current_header = ''
        current_sequence = ''
        # Switch telling if this the first entry
        virgin = True 
        for line in fasta_fh:
            # For usual headers
            if line[0] == '>' and virgin == False:
                yield(current_header, current_sequence)
                current_header = line[1:-1]
                current_sequence = ''
            # For the very first header of the file
            elif line[0] == '>' and virgin == True:
                virgin = False
                current_header = line[1:-1]
            # For sequence lines
            else:
                current_sequence += line[:-1]
        # For the last entry if file was not empty
        if not (current_header == '' and current_sequence == ''):
            yield(current_header, current_sequence)
