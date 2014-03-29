from subprocess import call
import os

class Lack(object):

    """A simple lack wrapper."""

    def __init__(self, lack_bin="lack", show_progress=False):
        self._lack_bin = lack_bin
        self._show_progress = show_progress

    def realign_reads(
        self, query_fasta_path, query_sam_path, ref_fasta_paths, 
            output_sam_path, nonmatch_path, accuracy=95, threads=1, 
            nonmatch_file=None, other_parameters=None):
        lack_call = [
            self._lack_bin, 
            "--query", query_sam_path,
            "--remapfilename", query_fasta_path,
            "--database"] + ref_fasta_paths + [
            "--outfile", output_sam_path,
            "--threads", str(threads),
            "--nomatchfilename", nonmatch_path]
        if self._show_progress is False:
            lack_call += ["--silent"]
        if other_parameters:
            lack_call.append(other_parameters)
        # Discard standard error output
        if self._show_progress is False:
            with open(os.devnull, "w") as devnull:
                call(lack_call, stderr=devnull)
        else:
            call(lack_call)            
