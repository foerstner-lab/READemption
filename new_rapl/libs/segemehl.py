from subprocess import call

class Segemehl(object):

    """A simple segemehl wrapper."""

    def __init__(self, segemehl_bin="segemehl"):
        self.segemehl_bin = segemehl_bin

    def build_index(self, fasta_files, index_file):
        """Create an index based on a list of fasta files"""
        call("%s -d %s -x %s" % 
             (self.segemehl_bin, " ".join(fasta_files), index_file),
             shell=True)

    def map_reads(read_file, index_file, output_file, nonmatch_file=None,
                  accurary=85.0, hit_strategy=1, other_parameters="",
                  sam_output=True, order=False, threads=1):
        pass
