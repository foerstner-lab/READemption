from subprocess import call

class Cutadapt(object):

    def __init__(self, args, cutadapt_options, cutadapt_bin="cutadapt"):
        self._cutadapt_bin = cutadapt_bin
        self._args = args 
        self._cutadapt_options = cutadapt_options
        
    def run_cutadapt_se(self, read_path, processed_read_path, lib_name):
        if self._cutadapt_options is None:
            cutadapt_options = []
            if self._args.adapter is not None:
                adapter_option = "-a " + self._args.adapter
                cutadapt_options.append(adapter_option)
            if self._args.min_read_length is not 12:
                min_read_length_option = "-m " + str(self._args.min_read_length)
                cutadapt_options.append(min_read_length_option)
            if self._args.min_phred_score is not None:
                min_phred_score_option = "-q " + str(self._args.min_phred_score)
                cutadapt_options.append(min_phred_score_option)
            if self._args.poly_a_clipping is True:
                poly_a_clipping_option = "-a AAAAAA" # Are 6 As a proper number
                                                     # for poly-A-clipping?
                cutadapt_options.append(poly_a_clipping_option)    
    
        else:
            cutadapt_options = self._cutadapt_options.split(",")
        output_path = "%s/%s_processed.fa.gz" % (
            processed_read_path, lib_name)
        cutadapt_call = [self._cutadapt_bin, "--quiet", "-o", output_path]
        for option in cutadapt_options:
            cutadapt_call.append(option)
        cutadapt_call.append(read_path)
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
