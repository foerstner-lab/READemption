import csv
import sys
import os
from subprocess import call

class DESeqRunner(object):

    def __init__(
            self, libs, conditions, deseq_raw_folder, deseq_extended_folder,
            deseq_script_path, gene_wise_quanti_combined_path, 
            deseq_tmp_session_info_script, deseq_session_info, 
            cooks_cutoff_off=False):
        self._libs = libs
        self._conditions = conditions
        self._deseq_raw_folder = deseq_raw_folder
        self._deseq_extended_folder = deseq_extended_folder
        self._deseq_script_path = deseq_script_path
        self._gene_wise_quanti_combined_path = gene_wise_quanti_combined_path
        self._deseq_tmp_session_info_script = deseq_tmp_session_info_script
        self._deseq_session_info = deseq_session_info
        self._cooks_cutoff_off = cooks_cutoff_off
        self._first_data_column = 11

    def write_session_info_file(self):
        with open(self._deseq_tmp_session_info_script, "w") as tmp_r_script_fh:
            tmp_r_script_fh.write("library('DESeq2')\nsessionInfo()\n")
        with open(self._deseq_session_info, "w") as session_info_fh:
            with open(os.devnull, "w") as devnull:
                call(["Rscript", self._deseq_tmp_session_info_script,], 
                     stdout=session_info_fh, stderr=devnull)
        os.remove(self._deseq_tmp_session_info_script)

    def create_deseq_script_file(self):
        libs_to_conditions = dict([
            (lib, condition) for lib, condition in
            zip(self._libs, self._conditions)])
        head_row = open(
            self._gene_wise_quanti_combined_path).readline()[:-1].split("\t")
        libs = head_row[self._first_data_column-1:]
        libs_str = ",".join(["'%s'" % lib for lib in libs])
        conditions = [libs_to_conditions[lib] for lib in libs]
        condition_str = ", ".join(["'%s'" % cond for cond in conditions])
        file_content = self._deseq_script_template() % (
            self._gene_wise_quanti_combined_path, self._first_data_column-1, 
            len(libs),self._first_data_column, libs_str, condition_str)
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
                if comparison_file_row[0] == "baseMean":
                    # Add another column to the header
                    comparison_file_row = [""] + comparison_file_row
                    # Extend column description
                    counting_file_row[self._first_data_column:] = [
                        "%s raw countings" % lib_name for lib_name in 
                        counting_file_row[self._first_data_column:]]
                output_fh.write("\t".join(
                    counting_file_row + comparison_file_row[1:]) + "\n")
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
        cooks_cutoff_str = ""
        if self._cooks_cutoff_off:
            cooks_cutoff_str = ", cooksCutoff=FALSE"
        for index, condition_combo in enumerate(condition_combos):
            call_string += ("comp%s <- results(dds, contrast="
                            "c('condition','%s', '%s')%s)\n" % (
                            index, condition_combo[0], condition_combo[1], 
                                cooks_cutoff_str))
            comparison_file = "deseq_comp_%s_vs_%s.csv" % (
                condition_combo[0], condition_combo[1])
            self._comparison_files_and_combos.append(
                (comparison_file, list(condition_combo)))
            call_string += (
                "write.table(comp%s, file='%s/%s', "
                "quote=FALSE, sep='\\t')\n" % (
                    index, self._deseq_raw_folder, comparison_file))
        return call_string

    def _deseq_script_template(self):
        return (
            "library('DESeq2')\n"
            "rawCountTable <- read.table('%s', skip=1, sep='\\t', "
            "quote='', comment.char='', colClasses=c(rep('character',%s), rep('numeric',%s)))\n"
            "countTable <- round(rawCountTable[,%s:length(names("
            "rawCountTable))])\n"
            "libs <- c(%s)\n"
            "conds <- c(%s)\n"
            "samples <- data.frame(row.names=libs, condition=conds)\n"
            "dds <- DESeqDataSetFromMatrix(countData=countTable, "
            "colData=samples, design=~condition)\n"
            "dds <- DESeq(dds)\n")
