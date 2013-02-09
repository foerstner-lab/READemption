from subprocess import call

class Segemehl(object):

    """A simple segemehl wrapper."""

    def __init__(self, segemehl_bin="segemehl"):
        self.segemehl_bin = segemehl_bin

    def build_index(self, fasta_files, index_file):
        """Create an index based on a list of fasta files"""
        call([self.segemehl_bin, "--database"] + fasta_files + [
             "--generate", index_file])

    def align_reads(self, read_file, index_file, fasta_files, output_file,
                  hit_strategy=1, accuracy=95, evalue=5.0, threads=1,
                  segemehl_format=False, order=False, nonmatch_file=None,
                  other_parameters=None):
        segemehl_call = [
            self.segemehl_bin, "--query", read_file,
            "--index", index_file,
            "--database"] + fasta_files + [
            "--outfile", output_file, 
            "--hitstrategy", str(hit_strategy),
            "--accuracy", str(accuracy),
            "--evalue", str(evalue),
            "--threads", str(threads)]
        if segemehl_format:
            segemehl_call.append("--SEGEMEHL")
        if order:
            segemehl_call.append("--order")
        if nonmatch_file:
            segemehl_call += ["--nomatchfilename", nonmatch_file]
        if other_parameters:
            segemehl_call.append(other_parameters)
        call(segemehl_call)

