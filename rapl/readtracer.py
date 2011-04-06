from rapl.fasta import FastaParser
from rapl.pathes import Pathes
from rapl.segemehl import SegemehlParser

class ReadTracer(object):

    def __init__(self):
        self.pathes = Pathes()
    
    def trace_reads(self):
        """Trace the way of each read during the different steps.

        A file is generated that can be used for downstream
        statistics.
        """
        for read_file in self.pathes.read_files:
            self.read_ids = []
            self.read_ids_and_traces = {}
            self._get_read_ids_and_lengths(read_file)
            self._read_first_mapping_output(read_file)
            self._read_first_mapping_unmapped_reads(read_file)
            self._read_clipped_unmapped_reads(read_file)
            self._read_size_filtered_reads_passed(read_file)
            self._read_size_filtered_reads_failed(read_file)
            self._read_second_mapping_output(read_file)
            self._read_second_mapping_unmapped_reads(read_file)
            self._read_combined_mapping_a_filtered_passed(read_file)
            self._read_combined_mapping_a_filtered_failed(read_file)
            self._write_trace_file(read_file)

    def _write_trace_file(self, read_file):
        """Write the trace of each read to a file.

        Arguments:
        - `read_file,`: the read file that was used to generate the
                        mapping files.

        """
        trace_fh = open(self.pathes.trace_file(read_file), "w")
        trace_fh.write("#Read id\tRead length\tNumber of mappings first run\t"
                       "length after clipping\tPassed size filter\t"
                       "Number of mappings second run\t"
                       "Passed a-content filter\tMapping length\t"
                       "Final status\n")
        for read_id in self.read_ids:
            trace = self.read_ids_and_traces[read_id]
            trace.setdefault("no_of_mappings_first_run", "-")
            trace.setdefault("length_after_clipping", "-")
            trace.setdefault("passed_size_filtering", "-")
            trace.setdefault("no_of_mappings_second_run", "-")
            trace.setdefault("passed_a-content_filtering", "-")
            trace.setdefault("mapping_length", "-")
            result_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                read_id, trace["length"], trace["no_of_mappings_first_run"],
                trace["length_after_clipping"], trace["passed_size_filtering"],
                trace["no_of_mappings_second_run"], 
                trace["passed_a-content_filtering"], 
                trace["mapping_length"],
                self._final_mapping_status(trace))
            trace_fh.write(result_line)

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
            self.pathes.read_file(read_file)):
            self.read_ids.append(header)
            self.read_ids_and_traces[header] = {'length' : len(seq)}

    def _read_first_mapping_output(self, read_file):
        """Read the result of the first Segemehl mapping.
        
        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the mapping file.

        """
        segemehl_parser = SegemehlParser()
        for entry in segemehl_parser.entries(
            self.pathes.raw_read_mapping_output(read_file)):
            entry_id = entry["id"][1:] # remove ">"
            self.read_ids_and_traces[entry_id].setdefault(
                "no_of_mappings_first_run", 0)
            self.read_ids_and_traces[entry_id]["no_of_mappings_first_run"] += 1

    def _read_first_mapping_unmapped_reads(self, read_file):
        """Read the file of unmapped reads of the first run.

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the mapping file.

        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self.pathes.unmapped_raw_read_file(read_file)):
                self.read_ids_and_traces[header]["no_of_mappings_first_run"] = 0

    def _read_clipped_unmapped_reads(self, read_file):
        """Read the file of clipped unmapped reads.

        Stores the length of the reads.

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the mapping file.

        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self.pathes.unmapped_read_clipped(read_file)):
            self.read_ids_and_traces[header][
                "length_after_clipping"] = len(seq)

    def _read_size_filtered_reads_passed(self, read_file):
        """Read the file of reads that passed the size filter.

        Arguments:
        - `read_file`: the name of the orignal read file

        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self.pathes.unmapped_clipped_size_filtered_read(read_file)):
            if header == "": continue
            self.read_ids_and_traces[header][
                "passed_size_filtering"] = True

    def _read_size_filtered_reads_failed(self, read_file):
        """Read the file of reads that did not pass the size filter.

        Arguments:
        - `read_file`: the name of the orignal read file

        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self.pathes.unmapped_clipped_size_failed_read(read_file)):
            if header == "": continue
            self.read_ids_and_traces[header][
                "passed_size_filtering"] = False

    def _read_second_mapping_output(self, read_file):
        """Read the file of the second Segemehl mapping.

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the mapping file.
        """
        segemehl_parser = SegemehlParser()
        for entry in segemehl_parser.entries(
            self.pathes.clipped_reads_mapping_output(read_file)):
            entry_id = entry["id"][1:] # remove ">"
            self.read_ids_and_traces[entry_id].setdefault(
                "no_of_mappings_second_run", 0)
            self.read_ids_and_traces[entry_id]["no_of_mappings_second_run"] += 1

    def _read_second_mapping_unmapped_reads(self, read_file):
        """Read the fasta file unmappable read of the second run.

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the first mapping file.
        """
        fasta_parser = FastaParser()
        for header, seq in fasta_parser.parse_fasta_file(
            self.pathes._unmapped_reads_second_mapping_path(read_file)):
            if header == "": continue
            self.read_ids_and_traces[header]["no_of_mappings_second_run"] = 0

    def _read_combined_mapping_a_filtered_passed(self, read_file):
        """Read file of mappings passing the A-contend filtering.

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the first mapping file.
        """
        segemehl_parser = SegemehlParser()
        for entry in segemehl_parser.entries(
            self.pathes.combined_mapping_file_a_filtered(read_file)):
            entry_id = entry["id"][1:] # remove ">"
            self.read_ids_and_traces[entry_id][
                "passed_a-content_filtering"] = True
            self.read_ids_and_traces[entry_id][
                "mapping_length"] = len(entry["sequence"])
    
    def _read_combined_mapping_a_filtered_failed(self, read_file):
        """Read file of mapping not passing the A-contend filtering.

        Arguments:
        - `read_file`: name of the read file that is used to generate
                       the first mapping file.
        """
        segemehl_parser = SegemehlParser()
        for entry in segemehl_parser.entries(
            self.pathes.combined_mapping_file_a_filter_failed(read_file)):
            entry_id = entry["id"][1:] # remove ">"
            self.read_ids_and_traces[entry_id][
                "passed_a-content_filtering"] = False

    def _final_mapping_status(self, trace):
        """Return the final mapping status of a read.

        Arguments:
        - `trace`: the trace of the a read.
        """
        if (trace["passed_a-content_filtering"] and 
            not trace["passed_a-content_filtering"] == "-"):
            if trace["no_of_mappings_first_run"] > 0:
                return("mapped_in_first_round")
            elif trace["no_of_mappings_second_run"] > 0:
                return("mapped_in_second_round")
        elif not trace["passed_a-content_filtering"]:
            if trace["no_of_mappings_first_run"] > 0:
                return("mapped_in_first_round-failed_a-content_filter")
            elif trace["no_of_mappings_second_run"] > 0:
                return("mapped_in_second_round-faild_a_content_filter")
        elif trace["passed_a-content_filtering"] == "-":
            if not trace["passed_size_filtering"]:
                return("failed_size_filter_after_clipping")
            if trace["no_of_mappings_second_run"] == 0:
                return("not_mappable_in_second_run")
        else:
            return("lost_somewhere")

    def create_tracing_summay(self):
        """
        """
        stati = [
            "mapped_in_first_round", 
            "mapped_in_first_round-failed_a-content_filter",
            "mapped_in_second_round", 
            "mapped_in_second_round-faild_a_content_filter",
            "failed_size_filter_after_clipping", "not_mappable_in_second_run",
            "lost_somewhere"]
        summary_fh = open(self.pathes.tracing_summary_file, "w")
        summary_fh.write(
            "#lib name\t" + 
            "\t".join(stati) +
            "\t" +
            "total number of reads\t" +
            "sum of mappable reads\t" + 
            "percentage mappable reads\t" + 
            "\n")
        for read_file in self.pathes.read_files:
            stati_and_countings = self._summarize_tracing_file(
                self.pathes.trace_file(read_file))
            countings = []
            for status in stati:
                stati_and_countings.setdefault(status, 0)
                countings.append(stati_and_countings[status])
            summary_fh.write(
                read_file + "\t" + 
                "\t".join([str(counting) for counting in countings]) + 
                "\t%s" % sum(countings) +
                "\t%s" % (stati_and_countings["mapped_in_first_round"] + 
                          stati_and_countings["mapped_in_second_round"]) +
                "\t%s" % ((stati_and_countings["mapped_in_first_round"] + 
                          stati_and_countings["mapped_in_second_round"])/
                          sum(countings)*100.0) +
                "\n")
        summary_fh.close()

    def _summarize_tracing_file(self, tracing_file):
        """
        """
        stati_and_countings = {}
        for line in open(tracing_file):
            if line[0] in ["#", "\n"]:
                continue
            final_status = line[:-1].split("\t")[8]
            stati_and_countings.setdefault(final_status, 0)
            stati_and_countings[final_status] += 1
        return(stati_and_countings)
