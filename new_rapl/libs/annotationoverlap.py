import sys
sys.path.append(".")
from libs.gff3 import Gff3Parser
from libs.intervaltree import IntervalTree, Interval
from libs.sam import SamParser

class AnnotationOverlap(object):

    """
    
    Find overlaps of genes and annotation independent of the strand
    and the length of the overlap. These raw overlaps will be filtered
    later.

    The following problems have to be solved: We have reads that are
    mapped to different reference genomes. For each reference genome
    there might be several annotation files. In principal an annotation
    files can have entries belonging to different reference
    genomes. The interval tree that is used to query for overlaps
    won't be able to handle annotation that have the same coordinates.

    """

    def __init__(self):
        # A dictionary of dictionaries of list
        # genomes -> intervals -> genes
        self.genome_interval_gene = {}
        self.genomes_and_intervals = {}
        self.genome_and_interval_trees = {}
        self.gff3_parser = Gff3Parser()
        self.sam_parser = SamParser()

    def read_annotation_files(self, annotation_file_paths):
        for annotation_file_path in annotation_file_paths:
            self._add_annotations(open(annotation_file_path))
        self._build_interval_trees()

    def search_overlaps(self, read_mapping_result_paths,
                        annotation_overlap_result_paths):
        for read_mapping_result_path, annotation_overlap_result_path  in zip(
            read_mapping_result_paths, annotation_overlap_result_paths):
            self._search_overlaps(
                open(read_mapping_result_path), 
                open(annotation_overlap_result_path, "w"))

    def _add_annotations(self, annotation_fh):
        for gff_entry in self.gff3_parser.entries(annotation_fh):
            self._add_anno_to_genome_interval_gene(gff_entry)
            self._add_anno_to_genome_and_interval(gff_entry)

    def _add_anno_to_genome_interval_gene(self, gff_entry):
        seq_id = gff_entry.seq_id
        interval_string = self._positions_to_key_string(
            gff_entry.start, gff_entry.end)
        self.genome_interval_gene.setdefault(seq_id, {})
        self.genome_interval_gene[seq_id].setdefault(interval_string, [])
        self.genome_interval_gene[seq_id][interval_string].append(gff_entry)
    
    def _add_anno_to_genome_and_interval(self, gff_entry):
        seq_id = gff_entry.seq_id
        self.genomes_and_intervals.setdefault(seq_id, [])
        self.genomes_and_intervals[
            seq_id].append((gff_entry.start, gff_entry.end))

    def _build_interval_trees(self):
        for seq_id, raw_intervals in self.genomes_and_intervals.items():
            intervals = []
            for raw_interval in raw_intervals:
                start, end = self._sorted_start_end(raw_interval[0], raw_interval[1])
                intervals.append(Interval(start, end))
            interval_tree = IntervalTree()
            interval_tree.fill(intervals)
            self.genome_and_interval_trees[seq_id] = interval_tree

    def _sorted_start_end(self, start, end):
        """Make sure that start and end are in the correct order."""
        return(sorted([start, end]))
            
    def _search_overlaps(self, sam_fh, output_fh):
        for sam_entry in self.sam_parser.entries(sam_fh):
            result_intervals = self._search_sam_entry_overlap(sam_entry)
            has_hit = False
            sam_seq_id = self._mod_seq_id(sam_entry.reference)
            used_intervals = {}
            for interval in result_intervals:
                interval_string = self._positions_to_key_string(
                    interval.start, interval.end)
                if interval_string in used_intervals:
                    continue
                used_intervals[interval_string] = 1
                for gene_gff in self._genes_of_interval(
                    sam_seq_id, interval.start, interval.end):
                    if sam_seq_id == gene_gff.seq_id:
                        has_hit = True
                        self._write_output_line_with_hit(
                            sam_entry, gene_gff, output_fh)
            if not has_hit:
                self._write_output_line_without_hit(sam_entry, output_fh)
            
    def _search_sam_entry_overlap(self, sam_entry):
        try:
            interval_tree = self.genome_and_interval_trees[sam_entry.reference]
        except KeyError:
            return([])
        start, end = self._sorted_start_end(sam_entry.start, sam_entry.end)
        interval = Interval(start, end)
        return(interval_tree.search_interval_overlaps(interval))

    def _genes_of_interval(self, seq_id, start, end):
        seq_id = self._mod_seq_id(seq_id)
        interval_string = self._positions_to_key_string(start, end)
        return(self.genome_interval_gene[seq_id][interval_string])

    def _write_output_line_with_hit(self, sam_entry, gene_gff, output_fh):
        output_fh.write(
            "\t".join([str(field) for field in
                       [sam_entry.query_id, sam_entry.reference,
                       sam_entry.start, sam_entry.end, sam_entry.strand,
                       sam_entry.number_of_hits_as_int, gene_gff.seq_id, 
                        gene_gff.source, gene_gff.feature, gene_gff.start, 
                        gene_gff.end, gene_gff.score, gene_gff.strand, 
                        gene_gff.phase, gene_gff.attribute_string]]) + "\n")

    def _write_output_line_without_hit(self, sam_entry, output_fh):
        output_fh.write(
            "\t".join([str(field) for field in 
                       [sam_entry.query_id, sam_entry.reference,
                       sam_entry.start, sam_entry.end, sam_entry.strand,
                       sam_entry.number_of_hits_as_int, "no_overlap"]]) + "\n")

    def _positions_to_key_string(self, start, end):
        return("%s-%s" % (start, end))

    def _mod_seq_id(self, seq_id):
        """ For files from NCBI. While the fasta file headers usually
        contain more information (e.g. "gi|15791399|ref|NC_002163.1|")
        the GFF files only use an accession number (in this case
        "NC_002163.1") for the refence genome.
        """
        if seq_id.startswith("gi|"):
            try:
                seq_id = seq_id.split("|")[3]
            except:
                pass
        return(seq_id)
