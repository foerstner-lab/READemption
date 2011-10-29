from subprocess import call

class Segemehl(object):

    """A simple segemehl wrapper."""

    def __init__(self, segemehl_bin="segemehl"):
        self.segemehl_bin = segemehl_bin

    def build_index(self, fasta_files, index_file):
        """Create an index based on a list of fasta files"""
        call("%s --database %s --generate %s " % 
             (self.segemehl_bin, " ".join(fasta_files), index_file),
             shell=True)

    def map_reads(self, read_file, index_file, fasta_files, output_file,
                  hit_strategy=1, accuracy=95, evalue=5.0, threads=1,
                  segemehl_format=False, order=False, nonmatch_file=None,
                  other_parameters=None):
        segemehl_call = (
            "%s --query %s --index %s --database %s --outfile %s " 
            "--hitstrategy %s --accuracy %s --evalue %s --threads %s"
            % (self.segemehl_bin, read_file, index_file, " ".join(fasta_files),
               output_file, hit_strategy, accuracy, evalue, threads))
        if segemehl_format:
            segemehl_call += " --SEGEMEHL"
        if order:
            segemehl_call += " --order"
        if nonmatch_file:
            segemehl_call += " --nomatchfilename %s" % nonmatch_file
        if other_parameters:
            segemehl_call += " " + other_parameters
        call(segemehl_call, shell=True)

