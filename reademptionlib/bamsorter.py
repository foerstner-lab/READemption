import pysam

class BamSorter(object):
    def __init__(
            self,
            tag,
            sort_by="name"
    ):
        self.tag = tag
        if sort_by == "name":
            self.sort_by = "-n"

    def sort_bam(self, input_bam, output_bam_sorted):
        pysam.sort(
            "-t",
            self.tag,
            self.sort_by,
            "-o",
            output_bam_sorted,
            input_bam,
        )