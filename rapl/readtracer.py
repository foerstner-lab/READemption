from libs.fasta import FastaParser
from libs.sam import SamParser
from rapl.parameters import Parameters
from rapl.paths import Paths

class ReadTracer(object):

    def __init__(self):
        self.paths = Paths()
        self.parameters = Parameters()        
    
    def trace_reads(self):
        """Trace the way of each read during the different steps.

        A file is generated that can be used for downstream
        statistics.
        """
        for read_file in self.paths.read_files:
            self.read_ids_and_traces = {}
            self._get_read_ids_and_lengths(read_file)
            self._read_clipped_reads_file(read_file)
            self._read_size_filtered_reads_passed_file(read_file)
            self._read_size_filtered_reads_failed_file(read_file)
            self._read_mapping_output_file(read_file)
            self._read_unmapped_reads_file(read_file)
            self._read_a_filtered_mapping_file(read_file)
            self._read_a_filter_failed_mapping_file(read_file)
            self._read_uniquely_mapped_reads_file(read_file)
            self._write_trace_file(read_file)

    def _write_trace_file(self, read_file):
        """Write the trace of each read to a file.

        Arguments:
        - `read_file,`: the read file that was used to generate the
                        mapping files.

        """
        trace_fh = open(self.paths.trace_file(read_file), "w")
        trace_fh.write("\t".join(
                ["#Read id", "Read length", "Number of mappings",
                 "Length after clipping", "Passed size filter", 
                 "Passed a-content filter", "Mapping length", 
                 "Passed Uniquely mapped filter", "Final status"]) + "\n")
        for read_id in self.read_ids_and_traces.keys():
            trace = self.read_ids_and_traces[read_id]
            trace.setdefault("no_of_mappings", "-")
            trace.setdefault("length_after_clipping", "-")
            trace.setdefault("passed_size_filtering", "-")
            trace.setdefault("passed_a-content_filtering", "-")
            trace.setdefault("mapping_length", "-")
            if self.parameters.uniquely_mapped_reads_only:
                trace.setdefault("passed_unique_mapping_filtering", False)
            else:
                trace.setdefault("passed_unique_mapping_filtering", "Not_checked")
            trace_fh.write("\t".join([str(x) for x in [
                        read_id, trace["length"], trace["no_of_mappings"],
                        trace["length_after_clipping"], trace["passed_size_filtering"],
                        trace["passed_a-content_filtering"], trace["mapping_length"],
                        trace["passed_unique_mapping_filtering"], 
                        self._final_mapping_status(trace)]]) + "\n")

    def _get_read_ids_and_lengths(self, read_file):
        """Get the ids and length of read sequences.

        Here all the initial reads are indexed and the length
        saved. The created dictionary is filled further in following
        steps.

        Arguments:
        - `read_file`: name of the read file
        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self.paths.read_file(read_file)):
            # TODO: TMP fix due to modification of the string
            # by segemehl
            header = self.mod_fasta_header(header)
            self.read_ids_and_traces[header] = {'length' : len(seq)}

    def _read_mapping_output_file(self, read_file):
        """Read the result of the first Sam mapping.
        
        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the mapping file.

        """
        sam_parser = SamParser()
        for entry in sam_parser.entries(
            self.paths.read_mapping_output_file(read_file)):
            self.read_ids_and_traces[entry["query"]].setdefault(
                "no_of_mappings", 0)
            self.read_ids_and_traces[entry["query"]]["no_of_mappings"] += 1

    def _read_unmapped_reads_file(self, read_file):
        """Read the file of unmapped reads of the first run.

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the mapping file.

        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self.paths.unmapped_reads_file(read_file)):
            # TODO: TMP fix due to modification of the string
            # by segemehl
            header = self.mod_fasta_header(header)
            self.read_ids_and_traces[header]["no_of_mappings"] = 0

    # TODO: TMP fix due to modification of the string by segemehl
    def mod_fasta_header(self, fasta_header):
        return(fasta_header.split("/")[0])

    def _read_clipped_reads_file(self, read_file):
        """Read the file of clipped unmapped reads.

        Stores the length of the reads.

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the mapping file.

        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self.paths.clipped_read_file(read_file)):
            # TODO: TMP fix due to modification of the string by segemehl
            header = self.mod_fasta_header(header)
            self.read_ids_and_traces[header][
                "length_after_clipping"] = len(seq)

    def _read_size_filtered_reads_passed_file(self, read_file):
        """Read the file of reads that passed the size filter.

        Arguments:
        - `read_file`: the name of the orignal read file

        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self.paths.clipped_size_filtered_read_file(
                read_file)):
            # TODO: TMP fix due to modification of the string by
            # segemehl
            header = self.mod_fasta_header(header)
            if header == "": continue
            self.read_ids_and_traces[header][
                "passed_size_filtering"] = True

    def _read_size_filtered_reads_failed_file(self, read_file):
        """Read the file of reads that did not pass the size filter.

        Arguments:
        - `read_file`: the name of the orignal read file

        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self.paths.clipped_size_filter_failed_read_file(read_file)):
            # TODO: TMP fix due to modification of the string by
            # segemehl
            header = self.mod_fasta_header(header)
            if header == "": continue
            self.read_ids_and_traces[header][
                "passed_size_filtering"] = False


    def _read_a_filtered_mapping_file(self, read_file):
        """Read file of mappings passing the A-contend filtering.

        Arguments:
        - `read_file`: name of the read file that is used to generat 
        the mapping file.
        """
        sam_parser = SamParser()
        for entry in sam_parser.entries(
            self.paths.a_filtered_mappings_file(read_file)):
            self.read_ids_and_traces[entry["query"]][
                "passed_a-content_filtering"] = True
            self.read_ids_and_traces[entry["query"]][
                "mapping_length"] = len(entry["sequence"])

            
    def _read_a_filter_failed_mapping_file(self, read_file):
        """Read file of mapping not passing the A-contend filtering.

        Arguments:
         `read_file`: name of the read file that is used to generate
                      the first mapping file.
        """
        sam_parser = SamParser()
        for entry in sam_parser.entries(
            self.paths.a_filter_failed_mappings_file(read_file)):
            self.read_ids_and_traces[entry["query"]][
                "passed_a-content_filtering"] = False


    def _read_uniquely_mapped_reads_file(self, read_file):
        if not self.parameters.uniquely_mapped_reads_only:
            return()
        sam_parser = SamParser()
        for entry in sam_parser.entries(
            self.paths.unique_mappings_only_file(read_file)):
            if self.parameters.uniquely_mapped_reads_only:
                self.read_ids_and_traces[entry["query"]][
                    "passed_unique_mapping_filtering"] = True

    def _final_mapping_status(self, trace):
        """Return the final mapping status of a read.

        Arguments:
        - `trace`: the trace of the a read.
        """
        if (trace["passed_a-content_filtering"] and 
            not trace["passed_a-content_filtering"] == "-"):
            if trace["no_of_mappings"] > 0:
                return("mapped")
        elif not trace["passed_a-content_filtering"]:
            if trace["no_of_mappings"] > 0:
                return("mapped_-failed_a-content_filter")
        elif trace["passed_a-content_filtering"] == "-":
            if not trace["passed_size_filtering"]:
                return("failed_size_filter_after_clipping")
        else:
            return("lost_somewhere")

    def create_tracing_summay(self):
        """
        """
        stati = ["mapped", "mapped_-failed_a-content_filter",
                 "failed_size_filter_after_clipping", "lost_somewhere"]
        summary_fh = open(self.paths.tracing_summary_file, "w")
        self._write_summary_header(stati, summary_fh)

        for read_file in self.paths.read_files:
            stati_and_countings, uniqely_mapped_read_countings = (
                self._summarize_tracing_file(self.paths.trace_file(read_file)))
            countings = []
            for status in stati:
                stati_and_countings.setdefault(status, 0)
                countings.append(stati_and_countings[status])
            summary_fh.write(
                "\t".join(
                    [read_file] + 
                    [str(number) for number in [
                            sum(countings), stati_and_countings["mapped"],
                            self._percentage_of_mapped_reads(
                                stati_and_countings, countings),
                            uniqely_mapped_read_countings,
                            self._percentage_of_uniquely_mapped_reads(
                                uniqely_mapped_read_countings, countings)]] +
                    [str(counting) for counting in countings]) + "\n")
        summary_fh.close()

    def _write_summary_header(self, stati, summary_fh):
        summary_fh.write("\t".join([
                    "#lib name", "total number of reads",  
                    "sum of mappable reads", "% mappable reads",
                    "uniquely mapped reads", "% of uniquely mapped reads"] 
                                   + stati) + "\n")

    def _percentage_of_mapped_reads(
        self, stati_and_countings, countings):
        try:
            return(
                round((stati_and_countings["mapped"]/sum(countings)*100.0),3))
        except ZeroDivisionError:
            return(0)

    def _percentage_of_uniquely_mapped_reads(
        self, uniqely_mapped_read_countings, countings):
        try:
            return(round((float(uniqely_mapped_read_countings) / 
                          sum(countings)*100.0),3))
        except ZeroDivisionError:
            return(0)

    def _summarize_tracing_file(self, tracing_file):
        """
        """
        stati_and_countings = {}
        uniquely_mapped_read_countings = 0
        for line in open(tracing_file):
            if line[0] in ["#", "\n"]:
                continue
            split_line = line[:-1].split("\t")
            final_status = split_line[8]
            no_of_mappings = split_line[2]
            if final_status == "mapped" and no_of_mappings == "1":
                uniquely_mapped_read_countings += 1
            stati_and_countings.setdefault(final_status, 0)
            stati_and_countings[final_status] += 1
        return(stati_and_countings, uniquely_mapped_read_countings)

