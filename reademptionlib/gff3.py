import csv
import sys


class Gff3Parser(object):
    """
    A format description can be found at:
    http://genome.ucsc.edu/FAQ/FAQformat.html#format3
    http://www.sequenceontology.org/gff3.shtml

    a validator can be found here:
    http://modencode.oicr.on.ca/cgi-bin/validate_gff3_online
    """

    def entries(self, input_gff_fh, annotation_name):
        """
        """
        gff_line_number = 0
        for entry_dict in csv.DictReader(
            input_gff_fh,
            delimiter="\t",
            fieldnames=[
                "seq_id",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ],
        ):
            gff_line_number += 1
            if entry_dict["seq_id"].startswith("#"):
                continue
            try:
                yield (self._dict_to_entry(entry_dict))
            except Exception as exception:
                sys.stderr.write(
                    f"Error in line number {gff_line_number} of {annotation_name}! Please make sure that you use GFF3 formated "
                    "annotation files. GTF2/GTF is not valid and the usage "
                    "of that format is not recommended anymore (see "
                    "http://www.sequenceontology.org/gff3.shtml for more "
                    "information).\n"
                    f"The line that caused the error was read in as: {entry_dict}"
                    f"The error message that occurred is: {exception} "
                )
                sys.exit(0)

    def _dict_to_entry(self, entry_dict):
        return Gff3Entry(entry_dict)


class Gff3Entry(object):
    def __init__(self, entry_dict):
        self.seq_id = entry_dict["seq_id"]
        self.source = entry_dict["source"]
        self.feature = entry_dict["feature"]
        # 1-based coordinates
        # Make sure that start <= end
        start, end = sorted([int(entry_dict["start"]), int(entry_dict["end"])])
        self.start = start
        self.end = end
        self.score = entry_dict["score"]
        self.strand = entry_dict["strand"]
        self.phase = entry_dict["phase"]
        self.attributes = self._attributes(entry_dict["attributes"])
        self.attribute_string = entry_dict["attributes"]

    def _attributes(self, attributes_string):
        """Translate the attribute string to dictionary"""
        if attributes_string is None:
            return {}
        if attributes_string.endswith(";"):
            attributes_string = attributes_string[:-1]
        return dict(
            [
                key_value_pair.split("=")
                for key_value_pair in attributes_string.split(";")
            ]
        )

    def __str__(self):
        return "\t".join(
            [
                str(field)
                for field in [
                    self.seq_id,
                    self.source,
                    self.feature,
                    self.start,
                    self.end,
                    self.score,
                    self.strand,
                    self.phase,
                    self.attribute_string,
                ]
            ]
        )
