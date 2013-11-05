from subprocess import call
import os

class Segemehl(object):

    """A simple segemehl wrapper."""

    def __init__(self, segemehl_bin="segemehl", show_progress=False):
        self._segemehl_bin = segemehl_bin
        self._show_progress = show_progress

    def build_index(self, fasta_files, index_file):
        """Create an index based on a list of fasta files"""
        segemehl_call = [self._segemehl_bin, "--database"] + fasta_files + [
            "--generate", index_file]
        if self._show_progress is False:
            with open(os.devnull, "w") as devnull:
                call(segemehl_call, stderr=devnull)
        else:
            call(segemehl_call)

    def align_reads(
        self, read_file_or_pair, index_file, fasta_files, output_file,
        hit_strategy=1, accuracy=95, evalue=5.0, threads=1, split=False, 
        segemehl_format=False, order=False, nonmatch_file=None,
        other_parameters=None, paired_end=False):
        if paired_end is False:
            assert type(read_file_or_pair) == str 
            segemehl_call = [
                self._segemehl_bin, 
                "--query", read_file_or_pair]
        else:
            assert type(read_file_or_pair) == list 
            segemehl_call = [
                self._segemehl_bin, 
                "--query", read_file_or_pair[0],
                "--mate", read_file_or_pair[1]]
        segemehl_call += [
            "--index", index_file,
            "--database"] + fasta_files + [
            "--outfile", output_file, 
            "--hitstrategy", str(hit_strategy),
            "--accuracy", str(accuracy),
            "--evalue", str(evalue),
            "--threads", str(threads)]
        if segemehl_format:
            segemehl_call.append("--SEGEMEHL")
        if order is True:
            segemehl_call.append("--order")
        if split is True:
            segemehl_call.append("--splits")
        if nonmatch_file:
            segemehl_call += ["--nomatchfilename", nonmatch_file]
        if self._show_progress is False:
            segemehl_call += ["--silent"]
        if other_parameters:
            segemehl_call.append(other_parameters)
        # Discard standard error output
        if self._show_progress is False:
            with open(os.devnull, "w") as devnull:
                call(segemehl_call, stderr=devnull)
        else:
            call(segemehl_call)
            
