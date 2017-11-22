from subprocess import call


class Cutadapt(object):
    def __init__(self, cutadapt_options, cutadapt_bin="cutadapt"):
        self._cutadapt_bin = cutadapt_bin
        self._cutadapt_options = cutadapt_options

    def run_cutadapt_se(self, read_path, processed_read_path, lib_name):
        output_path = "%s/%s_processed.fa.gz" % (
            processed_read_path, lib_name)
        cutadapt_call = [
            self._cutadapt_bin, "--quiet", "-o", output_path,
            self._cutadapt_options, read_path]
        call(cutadapt_call)

    def run_cutadapt_pe(
            self, read_path_pair, processed_read_path, lib_name):
        output_path_p1 = "%s/%s_p1_processed.fa.gz" % (
            processed_read_path, lib_name)
        output_path_p2 = "%s/%s_p2_processed.fa.gz" % (
            processed_read_path, lib_name)
        cutadapt_call = [self._cutadapt_bin, self._cutadapt_options,
                         "-o", output_path_p1, "-p",
                         output_path_p2, read_path_pair[0], read_path_pair[1]]
        call(cutadapt_call)
