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
        self._write_header()
        coutings = []
        for read_file in self.paths.read_files:
            self._count_and_write_results(read_file)

    def _count_and_write_results(self, read_file):
        no_of_mapped_reads, no_of_mappings, counting_by_fasta_header = (
            self._no_of_mappings_and_mapped_reads(read_file))
        # Get the countings per genome file via the fasta headers
        fasta_headers = [
            self.helper.get_header_of_genome_file(genome_file) 
            for genome_file in self.paths.genome_files]
        mapped_read_countings = [
            counting_by_fasta_header[fasta_header]["mapped_reads"] for
            fasta_header in fasta_headers]
        mapping_counting = [
            counting_by_fasta_header[fasta_header]["mappings"] for
            fasta_header in fasta_headers]
        combined_formated_countings = [
            str(round(counting, 3)) for counting in 
            [no_of_mapped_reads] + [no_of_mappings] +
            mapped_read_countings +  mapping_counting]
        self.summary_fh.write("\t".join(
                 [read_file] + combined_formated_countings) + "\n")

    def _write_header(self):
        self.summary_fh.write(
            "\t".join(
                ["#Library", "Total number of mapped reads", 
                 "Total number of mappings"] + 
                ["Mapped reads in " + genome_file 
                 for genome_file in self.paths.genome_files] +
                ["Mappings in " + genome_file 
                 for genome_file in self.paths.genome_files]) + "\n")

    def _no_of_mappings_and_mapped_reads(self, read_file):
        sam_parser = SamParser()
        fasta_headers_and_countings = {}
        # The total number of mapped reads could be calculated by
        # summing up the single countings. The counting of read by the
        # used method is a control that everything went fine.
        nonredundat_reads = {}
        number_of_mappings = 0
        # Initate dictionaries
        for genome_file in self.paths.genome_files:
            fasta_header = self.helper.get_header_of_genome_file(genome_file)
            fasta_headers_and_countings[fasta_header] = {
                "mapped_reads" : 0, "mappings" : 0}
        for entry in sam_parser.entries(
            self.paths.final_filtered_mapping_file(read_file)):
            fasta_headers_and_countings[
                entry["reference"]]["mapped_reads"] += (
                1.0/float(sam_parser.number_of_hits_as_int(entry)))
            fasta_headers_and_countings[
                entry["reference"]]["mappings"] += 1
            nonredundat_reads[entry["query"]] = 1
            number_of_mappings += 1 
        return(len(nonredundat_reads), number_of_mappings, 
               fasta_headers_and_countings)
