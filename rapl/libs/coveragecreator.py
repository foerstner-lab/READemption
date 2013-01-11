import pysam

class CoverageCreator(object):

    def __init__(self):
        self.elements_and_coverages = {"plus" : {}, "minus" : {}}

    def init_coverage_lists(self, bam_file):
        bam = pysam.Samfile(bam_file)
        for ref_seq, length in zip(bam.references, bam.lengths):
            for strand in ["plus", "minus"]:
                self.elements_and_coverages[strand][ref_seq] = [0.0] * length
        bam.close()

    def count_coverage(self, bam_file, read_count_splitting=True,
                       uniqueley_mapped_only=False):
        bam = pysam.Samfile(bam_file)
        for entry in bam.fetch():
            number_of_hits = dict(entry.tags)["NH"]
            if uniqueley_mapped_only and number_of_hits != 1:
                continue
            # Here a translation from 1-based system (SAM) to a
            # 0-based system (python lists) takes place. Due to this
            # each position is decreased by one. To cover the full
            # range of the end postion would need to be increased by
            # one. The substraction and addition result in a change of
            # zero.
            start = entry.pos - 1
            end = entry.aend
            ref_id = bam.getrname(entry.tid)

            # Normalize coverage increment by number of read mappings
            # per read
            if read_count_splitting:
                increment = 1.0 / float(number_of_hits)
            else:
                increment = 1.0
            if not entry.is_reverse:
                self.elements_and_coverages["plus"][ref_id][start:end] = [
                    coverage + increment for coverage in
                    self.elements_and_coverages["plus"][ref_id][start:end]]
            else:
                self.elements_and_coverages["minus"][ref_id][start:end] = [
                    coverage - increment for coverage in
                    self.elements_and_coverages["minus"][ref_id][start:end]]

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
