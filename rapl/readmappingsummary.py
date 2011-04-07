from libs.segemehl import SegemehlParser
from rapl.pathes import Pathes

class ReadMappingSummary(object):

    def __init__(self):
        self.pathes = Pathes()

    def create(self):
        """Create a file with lib/genome based read coutings."""
        coutings = []
        segemehl_parser = SegemehlParser()
        for read_file in self.pathes.read_files:
            counting_row = []
            for genome_file in self.pathes.genome_files:
                reads = {}
                for entry in segemehl_parser.entries(
                    self.pathes.combined_mapping_file_a_filtered_split(
                        read_file, genome_file)):
                    reads[entry["id"]] = 1
                counting_row.append(len(reads))
            coutings.append(counting_row)
        summary_fh = open(self.pathes.lib_genome_read_mapping_summary, "w")
        summary_fh.write("\t%s\n" % ("\t".join(self.pathes.genome_files)))
        for read_file, row in zip(self.pathes.read_files, coutings):
            summary_fh.write("%s\t%s\n" % (read_file, "\t".join(
                        [str(value) for value in row])))

