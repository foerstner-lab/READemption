from sam import SamParser, SamBuilder
import sys

def main():
    mapping_finder = MappingFinder()
    mapping_finder.read_overlap_files(sys.argv[3:])
    mapping_finder.read_mapping_file(sys.argv[1], sys.argv[2])

class MappingFinder(object):

    def read_overlap_files(self, overlap_files):
        self.mappings_with_overlaps = {}
        for overlap_file in overlap_files:
            self._read_overlap_file(open(overlap_file))

    def _read_overlap_file(self, overlap_fh):
        for line in overlap_fh:
            split_line = line.split("\t")
            mapping_key = "-".join(split_line[:3])
            self.mappings_with_overlaps[mapping_key] = 1

    def read_mapping_file(self, mapping_file, genome_accession):
        sam_parser = SamParser()
        sam_builder = SamBuilder()
        for entry in sam_parser.entries(mapping_file):
            if genome_accession != entry["reference"]
            start, end = sam_parser.entry_start_end_strand(entry)[:2]
            mapping_key = "-".join(
                [entry["query"], str(start), str(end)])
            if not mapping_key in self.mappings_with_overlaps:
                sys.stdout.write(sam_builder.entry_to_line(entry))

main()
