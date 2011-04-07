import os
from subprocess import PIPE
from subprocess import Popen
from libs.fasta import FastaParser
from rapl.pathes import Pathes

class InputStats(object):

    def __init__(self):
        self.pathes = Pathes()

    def create_read_file_stats(self):
        """Create a stat file for the input read files."""
        stat_fh = open(self.pathes.read_file_stats, "w")
        for read_file in self.pathes.read_files:
            stat_fh.write("%s:\n" % (read_file))
            stat_fh.write("* SHA256: %s\n" % (
                    self._sha256_of_file(self.pathes.read_file(read_file))))
            stat_fh.write("* Number of lines: %s\n" % (
                    self._number_of_lines_in_file(
                        self.pathes.read_file(read_file))))
            stat_fh.write("* Size: %s\n" % (
                    self._file_size(self.pathes.read_file(read_file))))
            stat_fh.write("* Number of Fasta entries: %s\n" % (
                    self._number_of_fasta_entries(
                        self.pathes.read_file(read_file))))
            stat_fh.write("\n")
        stat_fh.close()

    def create_genome_file_stats(self):
        """Create a stat file for the input genome files."""
        stat_fh = open(self.pathes.genome_file_stats, "w")
        for genome_file in self.pathes.genome_files:
            stat_fh.write("%s:\n" % (genome_file))
            stat_fh.write("* SHA256: %s\n" % (self._sha256_of_file(
                        self.pathes.genome_file(genome_file))))
            stat_fh.write("* Number of lines: %s\n" % (
                    self._number_of_lines_in_file(
                        self.pathes.genome_file(genome_file))))
            stat_fh.write("* Size: %s\n" % (self._file_size(
                                    self.pathes.genome_file(genome_file))))
            stat_fh.write("\n")
        stat_fh.close()

    def _sha256_of_file(self, file_path):
        """Calculate the SHA256 hash sum of a given file

        Arguments:
        - `file_path`: path of the file to process

        """
        # Todo: Fix handling of large files and then use this again
        # instead of calling the shell command
        #return(hashlib.sha256(open(file_path).read().encode()).hexdigest())
        return(Popen("sha256sum %s" % file_path , shell=True, stdout=PIPE
                     ).communicate()[0].split()[0].decode("utf-8"))

    def _number_of_lines_in_file(self, file_path):
        """Return the number of lines of a given file.

        Arguments:
        - `file_path`: path of the file to process

        """
        return(len(open(file_path).readlines()))

    def _file_size(self, file_path):
        """Return the size of a given file.

        Arguments:
        - `file_path`: path of the file to process

        """
        return(os.path.getsize(file_path))

    def _number_of_fasta_entries(self, file_path):
        """Return the number of fasta entries.

        Arguments:
        - `file_path`: path of the file to process
        
        """
        fasta_parser = FastaParser()
        counter = 0
        for head, seq in fasta_parser.parse_fasta_file(file_path):
            counter += 1
        return(counter)
