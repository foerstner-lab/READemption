import csv

class IntersectBedParser(object):

    def __init__(self, with_read_freq=False):
        self.with_read_freq = with_read_freq

    def entries(self, input_fh):
        fieldnames = ["sam_reference", "sam_mapping_start", "sam_mapping_end",
                      "sam_query_id", "sam_mapping_quality", "sam_strand",
                      "TODO1", "TODO2", "TODO3", "TODO4", "TODO5", "TODO6",
                      "gff_seq_id", "gff_source", "gff_feature", "gff_start",
                      "gff_end", "gff_score", "gff_strand", "gff_phase",
                      "gff_attribute_string", "overlap"]
        if self.with_read_freq:
            fieldnames.appedend("read_mapping_freq")
        for entry_dict in csv.DictReader(input_fh, delimiter="\t",
            fieldnames=fieldnames):
            yield(self._dict_to_entry(entry_dict))

    def _dict_to_entry(self, entry_dict):
        return(IntersectionEntry(entry_dict))

class IntersectionEntry(object):

    def __init__(self, entry_dict):
        self.sam_reference = entry_dict["sam_reference"]
        self.sam_start = int(entry_dict["sam_mapping_start"])
        self.sam_end = int(entry_dict["sam_mapping_end"])
        self.sam_query_id = entry_dict["sam_query_id"]
        self.sam_mapping_quality = entry_dict["sam_mapping_quality"]
        self.sam_strand = entry_dict["sam_strand"]
        self.gff_seq_id = entry_dict["gff_seq_id"]
        self.gff_source = entry_dict["gff_source"]
        self.gff_feature = entry_dict["gff_feature"]
        self.gff_start = int(entry_dict["gff_start"])
        self.gff_end = int(entry_dict["gff_end"])
        self.gff_score = entry_dict["gff_score"]
        self.gff_strand = entry_dict["gff_strand"]
        self.gff_phase = entry_dict["gff_phase"]
        self.gff_attribute_string = entry_dict["gff_attribute_string"]
        self.overlap = int(entry_dict["overlap"])
        self.orientation = self._orientation()
        if "read_mapping_freq" in entry_dict:
            self.read_mapping_freq = int(entry_dict["read_mapping_freq"])

    def _orientation(self):
        if self.sam_strand == self.gff_strand:
            return("s") # sense
        else:
            return("a") # anti-sense
