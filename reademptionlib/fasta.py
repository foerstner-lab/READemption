class FastaParser(object):
    """ A class for parsing and handling fasta files. """

    def entries(self, fasta_fh):
        """Go line by line though the file and return fasta entries. """
        current_header = ""
        current_sequence = ""
        # Switch telling if this the first entry
        virgin = True
        for line in fasta_fh:
            # For usual headers
            if line[0] == ">" and virgin is False:
                yield (current_header, current_sequence)
                current_header = line[1:-1]
                current_sequence = ""
            # For the very first header of the file
            elif line[0] == ">" and virgin is True:
                virgin = False
                current_header = line[1:-1]
            # For sequence lines
            else:
                current_sequence += line[:-1]
        # For the last entry if file was not empty
        if not (current_header == "" and current_sequence == ""):
            yield (current_header, current_sequence)

    def single_entry_file_header(self, fasta_fh):
        first_line = fasta_fh.readline()
        header = first_line[1:-1]
        fasta_fh.seek(0)
        return header

    def header_id(self, header):
        """Return only the id of a fasta header

        Everything after the first white space is discard.
        """
        return header.split()[0]
