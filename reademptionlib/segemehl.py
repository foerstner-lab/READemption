from subprocess import call
import os


class Segemehl(object):

    def __init__(self, segemehl_bin="segemehl", show_progress=False):
        self._segemehl_bin = segemehl_bin
        self._show_progress = show_progress

    def build_index(self, fasta_files, index_file):
        segemehl_call = [
            self._segemehl_bin, "--database"] + fasta_files + [
            "--generate", index_file]
        if self._show_progress is False:
            with open(os.devnull, "w") as devnull:
                call(segemehl_call, stderr=devnull)
        else:
            call(segemehl_call)

    def run_alignment(
        self, read_file_or_pair, index_file, fasta_files, output_file,
            threads, nomatch_path, hit_strategy, accuracy=95, evalue=5.0,
            split=False, segemehl_format=False, order=False,
            other_parameters=None, paired_end=False):
        if not paired_end:
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
            "--threads", str(threads), 
            "--nomatchfilename", nomatch_path,
            "--hitstrategy", str(hit_strategy),
            "--accuracy", str(accuracy),
            "--evalue", str(evalue)]
        if order is True:
            segemehl_call.append("--order")
        if split is True:
            segemehl_call.append("--splits")
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
