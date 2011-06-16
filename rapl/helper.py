from rapl.paths import Paths

class Helper(object):
    """Offers misc help for other classes"""

    def __init__(self):
        self.paths = Paths()

    def get_header_of_genome_file(self, genome_file):
        """Return the shor header of a given fasta file."""
        genome_file_path = self.paths.genome_file(genome_file)
        return(open(genome_file_path).readline().split()[0][1:])

