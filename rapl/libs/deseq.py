import csv
import sys
from subprocess import call

class DESeqRunner(object):

    def __init__(
            self, libs, conditions, deseq_raw_folder, deseq_extended_folder,
            deseq_script_path, gene_wise_quanti_combined_path, 
            no_replicates=False):
        self._libs = libs
        self._conditions = conditions
        self._deseq_raw_folder = deseq_raw_folder
        self._deseq_extended_folder = deseq_extended_folder
        self._deseq_script_path = deseq_script_path
        self._gene_wise_quanti_combined_path = gene_wise_quanti_combined_path
        self._no_replicastes = no_replicates
        self._first_data_column = 11

    def create_deseq_script_file(self):
        libs_to_conditions = dict([
            (lib, condition) for lib, condition in
            zip(self._libs, self._conditions)])
        head_row = open(
            self._gene_wise_quanti_combined_path).readline()[:-1].split("\t")
        libs = head_row[self._first_data_column-1:]
        conditions = [libs_to_conditions[lib] for lib in libs]
        condition_str = ", ".join(["'%s'" % cond for cond in conditions])
        dispersion_params = ""
        if self._no_replicastes is True:
            dispersion_params = ", method='blind', sharingMode='fit-only'"
        file_content = self._deseq_script_template() % (
            self._gene_wise_quanti_combined_path, self._first_data_column,
            condition_str, dispersion_params)
        file_content += self._comparison_call_strings(conditions)
        deseq_fh = open(self._deseq_script_path, "w")
        deseq_fh.write(file_content)
        deseq_fh.close()

    def run_deseq(self):
        call(["Rscript", self._deseq_script_path])

    def merge_counting_files_with_results(self):
        for comparison_file, combo in self._comparison_files_and_combos:
            output_fh = open("%s/%s" % (
                self._deseq_extended_folder,
                comparison_file.replace(
                ".csv", "_with_annotation_and_countings.csv")), "w")
            try:
                deseq_result_fh = open("%s/%s" % (
                    self._deseq_raw_folder, comparison_file))
            except:
                sys.stderr.write("Apparently DESeq did not generate the "
                                 "file \"%s\". Extension stopped.\n" % 
                                 comparison_file)
                continue
            for counting_file_row, comparison_file_row in zip(
                csv.reader(open(
                    self._gene_wise_quanti_combined_path), delimiter="\t"),
                csv.reader(deseq_result_fh, delimiter="\t")):
                if comparison_file_row[0] == "id":
                    comparison_file_row = [""] + comparison_file_row
                    # Add condition name to the headline
                    comparison_file_row[3] += " (%s)" % combo[0]
                    comparison_file_row[4] += " (%s)" % combo[1]
                output_fh.write("\t".join(
                    counting_file_row + comparison_file_row[2:]) + "\n")
            output_fh.close()

    def _condition_combos(self, conditions):
        for cond1 in conditions[:]:
            for cond2 in conditions[:]:
                if not cond1 == cond2:
                    yield((cond1, cond2))

    def _comparison_call_strings(self, conditions):
        call_string = ""
        condition_combos = self._condition_combos(conditions)
        self._comparison_files_and_combos = []
        for index, condition_combo in enumerate(condition_combos):
            call_string += "comp%s <- nbinomTest(cds, '%s', '%s')\n" % (
                index, condition_combo[0], condition_combo[1])
            comparison_file = "deseq_comp_%s_vs_%s.csv" % (
                condition_combo[0], condition_combo[1])
            self._comparison_files_and_combos.append(
                (comparison_file, list(condition_combo)))
            call_string += (
                "write.table(comp%s, file='%s/%s', "
                "quote=FALSE, sep='\\t')\n" % (
                    index, self._deseq_raw_folder, comparison_file))
        return(call_string)

    def _deseq_script_template(self):
        return(
            "library('DESeq')\n"
            "rawCountTable <- read.table('%s', skip=1, sep='\\t')\n"
            "names(rawCountTable)\n"
            "countTable <- round(rawCountTable[,%s:length(names(rawCountTable))])\n"
            "conds <- c(%s)\n"
            "cds <- newCountDataSet(countTable, conds)\n"
            "cds <- estimateSizeFactors(cds)\n"
            "cds <- estimateDispersions(cds%s)\n")
