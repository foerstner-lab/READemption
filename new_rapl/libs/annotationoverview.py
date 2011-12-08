import csv

class AnnotationOverview(object):

    # TODO 
    # - Read the overlap tables 
    #   - 1 count the number of overlaps of each mapping with how many gene
    #   - 2 use these countings in a second run
    # - Once all counting are done read each GFF file and print for each line the 
    #   coutings

    def __init__(self, needed_overlap=1):
        
        # This dictionary contain the final results. The countings are stored by
        # read file, orientation (sense/antisense) and gene key.
        # self.annotation_hit_countings[read_file][orientation][gene_key]
        self.annotation_hit_countings = {}
        self.mappings_and_no_of_overlaps = {}
        self.annotation_overlap_parser = AnnotationOverlapParser()

    def read_annotation_overlap_files(self):
        for annotation_overlap_file, read_file_name in zip():
            self._read_annotation_overlap_file(
                open(annotation_overlap_file), read_file_name)

    def _read_annotation_overlap_file(
        self, annotation_overlap_fh, read_file_name):
        """
        The annotation file is read twice. The first time the number
        of valid overlaps for each read mapping is counted. For this
        both directions - sense and antisense - are taken into
        account. The second time these numbers of total overlap of each
        mapping are used to normalize the each overlap counting.
        
        The sense and antisense values are currently counted
        separately and will be add together. At late developing stages
        this could be used to modify the normalization if wanted.
        """
        for annotation_overlap in self.annotation_overlap_parser.entries(
            annotation_overlap_fh):
            self._count_overlaps_per_mapping(annotation_overlap)
            
        for annotation_overlap in self.annotation_overlap_parser.parse(
            annotation_overlap_fh):
            pass

        # The values were only valid for the current file overlapping
        # file and only useful until now. So the dictionary needs to
        # be reset now.
        self.mappings_and_no_of_overlaps.clear()
            
    def _count_overlaps(self, annotation_overlap, read_file_name):
        if annotation_overlap.overlap_length < self.needed_overlap:
            return
        gene_key_string = self._gene_key_string(annotation_overlap)
        # Take care that all dictionaries have required keys
        self.annotation_hit_countings.setdefault(read_file_name, {})
        self.annotation_hit_countings[read_file_name].setdefault(
            annotation_overlap.orientation, {})
        self.annotation_hit_countings[read_file_name][
            annotation_overlap.orientation].setdefault(gene_key_string, {})
        # Add the value
        self.annotation_hit_countings[read_file_name][
            annotation_overlap.orientation][
            gene_key_string] += self._value_to_add(annotation_overlap)
        
    def _value_to_add(self, annotation_overlap):
        mapping_key_string = self._mapping_key_string(annotation_overlap)
        number_of_mappings = annotation_overlap.sam_number_of_hits
        # Sum up the number of overlapping genes for sense and antisense
        number_of_gene_overlaps = sum(self.mappings_and_no_of_overlaps[
            mapping_key_string].items())
        return(1.0 / number_of_mappings / number_of_gene_overlaps )

    def _count_overlaps_per_mapping(self, annotation_overlap):
        if annotation_overlap.overlap_length < self.needed_overlap:
            return
        mapping_key_string = self._overlap_key_string(annotation_overlap)
        self.mappings_and_no_of_overlaps.setdefault(mapping_key_string, {})
        self.mappings_and_no_of_overlaps[
            mapping_key_string][annotation_overlap.orientation] += 1

    def _mapping_key_string(self, annotation_overlap):
        pass
    
    def _gene_key_string(self, annotation_overlap):
        pass
    
    def print_result_table(self, annotation_file_paths, read_file_names):
        for annotation_file_path in annotation_file_paths:
            self._print_result_table(annotation_file_path, read_file_names)

    def _print_result_table(self, annotation_file_path, read_file_names, output_fh):
        gff3_parser = Gff3Parser()
        self._write_header()
        for gff3_entry in gff3_parser.entries(open(annotation_file_path)):
            gene_key_string = self._gff3_entry_gene_key_string(gff3_entry)
            values = [self.annotation_hit_countings[read_file][orientation].get(
                    gene_key_string, 0) 
                      for read_file in read_file_names 
                      for orientation in ["s", "a"]]
            output_fh.write(gff3_entry.string() + "\t" + "\t".join(values))

    def _gff3_entry_gene_key_string(self, gff3_entry):
        pass
        
    def _write_header(self):
        pass

class AnnotationOverlapParser(object):

    def entries(self, input_fh):
        for row in csv.DictReader(input_fh, delimiter="\t", fieldnames=[
                "sam_query_id", "sam_reference", "sam_start", "sam_end", 
                "sam_strand", "sam_number_of_hits", "gff_reference", 
                "gff_source", "gff_feature", "gff_start", "gff_end", 
                "gff_score", "gff_strand", "ggf_phase", "gff_attibutes"]):
            yield(self._dict_to_entry(row))

    def _dict_to_entry(self, entry_dict):
        return(AnnotationOverlap(entry_dict))

class AnnotationOverlap(object):

    def __init__(self, entry_dict):
        self.sam_query_id = entry_dict["sam_query_id"]
        self.sam_reference = entry_dict["sam_reference"]
        self.sam_start = int(entry_dict["sam_start"])
        self.sam_end = int(entry_dict["sam_end"])
        self.sam_strand = entry_dict["sam_strand"]
        self.sam_number_of_hits = int(entry_dict["sam_number_of_hits"])
        self.gff_reference = entry_dict["gff_reference"]
        self.gff_source = entry_dict["gff_source"]
        self.gff_feature = entry_dict["gff_feature"]
        self.gff_start = int(entry_dict["gff_start"])
        self.gff_end = int(entry_dict["gff_end"])
        self.gff_score = entry_dict["gff_score"]
        self.gff_strand = entry_dict["gff_strand"]
        self.ggf_phase = entry_dict["ggf_phase"]
        self.gff_attibutes = entry_dict["gff_attibutes"]
        self.orientation = self._orientation()
        self.overlap_length = self._overlap_length()
        
    def _orientation(self):
        if self.sam_strand == self.gff_strand:
            return("s") # sense
        else:
            return("a") # anti-sense

    def _overlap_length(self):
        return(min([self.sam_end, self.gff_end]) - 
               max([self.sam_start, self.gff_start]) + 1)
