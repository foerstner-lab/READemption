from rapl.paths import Paths

class ReadTracerViz(object):
    
    def __init__(self):
        self.paths = Paths()

    def create_mapping_length_histograms(self):
        r_file = open(self.paths.mapping_length_hist_r_file, "w")
        r_file.write("library(\"ggplot2\")\n")
        r_file.write("pdf(\"../../%s\")\n" % self.paths.mapping_length_hist_pdf_file)
        for read_file, index in zip(
            self.paths.read_files, range(len(self.paths.read_files))):
            r_file.write(self._create_mapping_length_hist_string(read_file, index))

    def _create_mapping_length_hist_string(self, read_file, index):
        return("trace_data_%s <- read.table(\"../../%s\", skip=1, "
               "na.string=\"-\", comment.char=\"\")\n"
               "qplot(trace_data_%s$V7, data=trace_data_%s, binwidth=1, "
               "xlab=\"Mapping length [nt]\", "
               "main=\"Mapping length distribution for %s\")\n" % (
                index, self.paths.trace_file(read_file), index, index, read_file))        
