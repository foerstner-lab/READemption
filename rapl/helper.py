class Helper(object):
    """Offers misc help for other classes"""

    def get_header_from_fasta_file(self, fasta_file):
        """Return the shor header of a given fasta file."""
        return(open(fasta_file).readline().split()[0][1:])

