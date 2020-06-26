from reademptionlib.segemehl import Segemehl


class ReadAligner(object):
    """An abstraction layer for different short read aligners."""

    def __init__(self, segemehl_bin, show_progress):
        self.segemehl = Segemehl(segemehl_bin, show_progress=show_progress)

    def build_index(self, ref_seq_paths, index_path):
        self.segemehl.build_index(ref_seq_paths, index_path)

    def run_alignment(
        self,
        read_path_or_pair,
        index_path,
        ref_seq_path,
        output_path,
        nomatch_path,
        threads,
        accuracy,
        evalue,
        split,
        paired_end=False,
    ):
        self.segemehl.align_reads(
            read_path_or_pair,
            index_path,
            ref_seq_path,
            output_path,
            nonmatch_file=nomatch_path,
            threads=threads,
            accuracy=accuracy,
            evalue=evalue,
            split=split,
            paired_end=paired_end,
        )
