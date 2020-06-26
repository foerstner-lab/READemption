class FastqParser(object):
    """A class for parsing and handling fastQ files.

    Currently only taking care of the sequences not the quality
    strings. Only four line entries can be handled
    """

    def entries(self, fastq_fh):
        """Go line by line though the file and return fastq entries. """
        current_header = ""
        current_sequence = ""
        # Switch telling if this the first entry
        virgin = True
        entry_line_counter = 0
        for line in fastq_fh:
            # For usual headers
            if line[0] == "@" and virgin is False and entry_line_counter == 4:
                entry_line_counter = 1
                yield (current_header, current_sequence)
                current_header = line[1:-1]
                current_sequence = ""
            # For the very first header of the file
            elif line[0] == "@" and virgin is True:
                virgin = False
                current_header = line[1:-1]
                entry_line_counter = 1
            # For sequence lines
            else:
                entry_line_counter += 1
                if entry_line_counter == 2:
                    current_sequence = line[:-1]
        # For the last entry if file was not empty
        if not (current_header == "" and current_sequence == ""):
            yield (current_header, current_sequence)

    def single_entry_file_header(self, fastq_fh):
        first_line = fastq_fh.readline()
        header = first_line[1:-1]
        fastq_fh.seek(0)
        return header

    def header_id(self, header):
        """Return only the id of a fastq header

        Everything after the first white space is discard.
        """
        return header.split()[0]
