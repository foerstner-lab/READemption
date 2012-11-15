from libs.sam import SamParser

class CoverageCreator(object):

    def __init__(self, samtools_bin="samtools"):
        self._sam_parser = SamParser(samtools_bin)
        self.elements_and_coverages = {"plus" : {}, "minus" : {}}

    def init_coverage_lists(self, bam_file):
        for ref_seq, length in self._sam_parser.ref_seq_ids_and_lengths_bam(
            bam_file).items():
            for strand in ["plus", "minus"]:
                self.elements_and_coverages[strand][ref_seq] = [0.0] * length

    def count_coverage(self, bam_file, read_count_splitting=True,
                       uniqueley_mapped_only=False):
        for entry in self._sam_parser.entries_bam(bam_file):
            if uniqueley_mapped_only and entry.number_of_hits_as_int != 1:
                continue
            # Here a translation from 1-based system (SAM) to a
            # 0-based system (python lists) takes place. Due to this
            # each position is decreased by one. To cover the full
            # range of the end postion would need to be increased by
            # one. The substraction and addition result in a change of
            # zero.
            start = entry.start - 1
            end = entry.end
            # Normalize coverage increment by number of read mappings
            # per read
            if read_count_splitting:
                increment = 1.0 / float(entry.number_of_hits_as_int)
            else:
                increment = 1.0
            if entry.strand == "+":
                self.elements_and_coverages["plus"][entry.reference][
                    start:end] = [
                    coverage + increment for coverage in
                    self.elements_and_coverages["plus"][entry.reference][
                            start:end]]
            elif entry.strand == "-":
                self.elements_and_coverages["minus"][entry.reference][
                    start:end] = [
                    coverage - increment for coverage in
                    self.elements_and_coverages["minus"][entry.reference][
                            start:end]]

    def write_to_files(self, output_file_prefix, read_file_name, factor=1.0,
                      output_format="wiggle"):
        if output_format == "wiggle":
            self._write_to_wiggle_files(
                output_file_prefix, read_file_name, factor=factor)

    def _write_to_wiggle_files(self, output_file_prefix, read_file_name, factor):
        for strand in ["plus", "minus"]:
            output_fh = open("%s_%s.wig" % (output_file_prefix, strand), "w")
            output_fh.write("track type=wiggle_0 name=\"%s_%s\"\n" % (
                    read_file_name, strand))
            for element in sorted(self.elements_and_coverages[strand].keys()):
                output_fh.write("variableStep chrom=%s span=1\n" % (element))
                # Filter values of 0 and multiply other the remaining
                # ones by the given factor. pos is increased by 1 as a
                # translation from a 0-based sysem (Python list) to a
                # 1 based system (wiggle) takes place.
                output_fh.write(
                    "\n".join(
                        ["%s %s" % (pos + 1, coverage * factor)
                         for pos, coverage in
                         filter(lambda pos_and_cov: pos_and_cov[1] != 0.0,
                                enumerate(self.elements_and_coverages[
                                        strand][element]))]) + "\n")
            output_fh.close()

