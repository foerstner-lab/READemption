from libs.sam import SamParser
from rapl.paths import Paths

class ReadMappingSummary(object):

    def __init__(self):
        self.paths = Paths()

    def create(self):
        """Create a file with lib/genome based read coutings."""
        coutings = []
        sam_parser = SamParser()
        for read_file in self.paths.read_files:
            counting_row = []
            for genome_file in self.paths.genome_files:
                reads = {}
                for entry in sam_parser.entries(
                    self.paths.combined_mapping_file_a_filtered_split(
                        read_file, genome_file)):
                    reads[entry["query"]] = 1
                counting_row.append(len(reads))
            coutings.append(counting_row)
        summary_fh = open(self.paths.lib_genome_read_mapping_summary, "w")
        summary_fh.write("\ttotal\t%s\n" % ("\t".join(self.paths.genome_files)))
        for read_file, row in zip(self.paths.read_files, coutings):
            summary_fh.write("%s\t%s\t%s\n" % (
                    read_file,
                    sum(row),
                    "\t".join([str(value) for value in row])))
