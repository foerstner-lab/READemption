from libs.sam import SamParser
from rapl.paths import Paths
from rapl.helper import Helper

class ReadMappingSummary(object):

    def __init__(self):
        self.paths = Paths()
        self.helper = Helper()
        self.summary_fh = open(self.paths.lib_genome_read_mapping_summary, "w")

    def create(self):
        """Create a file with lib/genome based read coutings."""

        coutings = []
        self.summary_fh.write("\ttotal\t%s\n" % ("\t".join(self.paths.genome_files)))
        for read_file in self.paths.read_files:
            total_counting, counting_by_fasta_header = (
                self._no_of_reads_mapped_to_genome_file(read_file))
            # Get the number of mapped read of per genome file via the
            # fasta headers
            countings = [
                counting_by_fasta_header.get(
                    self.helper.get_header_of_genome_file(genome_file), 0)
                for genome_file in self.paths.genome_files]
            # Format
            countings = [str(round(counting, 3)) for counting in countings]
            self.summary_fh.write("\t".join(
                    [read_file] + [str(total_counting)] + countings) + "\n")
            
    def _no_of_reads_mapped_to_genome_file(self, read_file):
        sam_parser = SamParser()
        fasta_headers_and_countings = {}
        # The total number of mapped reads could be calculated by
        # summing up the single countings. The counting of read by the
        # used method is a control that everything went fine.
        nonredundat_reads = {}
        for entry in sam_parser.entries(
            self.paths.combined_mapping_file_a_filtered(read_file)):
            fasta_headers_and_countings.setdefault(entry["reference"], 0)
            fasta_headers_and_countings[entry["reference"]] += (
                1.0/float(sam_parser.number_of_hits_as_int(entry)))
            nonredundat_reads[entry["query"]] = 1
        return(len(nonredundat_reads), fasta_headers_and_countings)
