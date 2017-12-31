import csv
import sys
import os
import pandas as pd
from subprocess import call


class DESeq2Runner(object):

    def __init__(
            self, libs, conditions, deseq2_raw_folder, deseq2_extended_folder,
            deseq2_script_path, deseq2_pca_heatmap_path,
            gene_wise_quanti_combined_path,
            deseq2_tmp_session_info_script, deseq2_session_info,
            cooks_cutoff_off=False):
        self._libs = libs
        self._conditions = conditions
        self._deseq2_raw_folder = deseq2_raw_folder
        self._deseq2_extended_folder = deseq2_extended_folder
        self._deseq2_script_path = deseq2_script_path
        self._deseq2_pca_heatmap_path = deseq2_pca_heatmap_path
        self._gene_wise_quanti_combined_path = gene_wise_quanti_combined_path
        self._deseq2_tmp_session_info_script = deseq2_tmp_session_info_script
        self._deseq2_session_info = deseq2_session_info
        self._cooks_cutoff_off = cooks_cutoff_off
        self._first_data_column = 11

    def write_session_info_file(self):
        with open(self._deseq2_tmp_session_info_script, "w") as tmp_r_script_fh:
            tmp_r_script_fh.write("library('DESeq2')\nsessionInfo()\n")
        with open(self._deseq2_session_info, "w") as session_info_fh:
            with open(os.devnull, "w") as devnull:
                call(["Rscript", self._deseq2_tmp_session_info_script],
                     stdout=session_info_fh, stderr=devnull)
        os.remove(self._deseq2_tmp_session_info_script)

    def create_deseq2_script_file(self):
        libs_to_conditions = dict([
            (lib, condition) for lib, condition in
            zip(self._libs, self._conditions)])
        head_row = open(
            self._gene_wise_quanti_combined_path).readline()[:-1].split("\t")
        libs = head_row[self._first_data_column-1:]
        libs_str = ",".join(["'%s'" % lib for lib in libs])
        conditions = [libs_to_conditions[lib] for lib in libs]
        condition_str = ", ".join(["'%s'" % cond for cond in conditions])
        file_content = self._deseq2_script_template() % (
            self._gene_wise_quanti_combined_path, self._first_data_column-1,
            len(libs), self._first_data_column, libs_str, condition_str,
            self._deseq2_pca_heatmap_path)
        file_content += self._comparison_call_strings(conditions)
        deseq2_fh = open(self._deseq2_script_path, "w")
        deseq2_fh.write(file_content)
        deseq2_fh.close()

    def run_deseq2(self):
        call(["Rscript", self._deseq2_script_path])

    def merge_counting_files_with_results(self):
        for comparison_file, combo in self._comparison_files_and_combos:
            output_fh = open("%s/%s" % (
                self._deseq2_extended_folder,
                comparison_file.replace(
                    ".csv", "_with_annotation_and_countings.csv")), "w")
            output_fh.write("# Reference library (divisor): {}\n".format(
                combo[1]))
            output_fh.write("# Comparison library (numerator): {}\n".format(
                combo[0]))
            try:
                deseq2_result_fh = open("%s/%s" % (
                    self._deseq2_raw_folder, comparison_file))
            except:
                sys.stderr.write("Apparently DESeq2 did not generate the "
                                 "file \"%s\". Extension stopped.\n" %
                                 comparison_file)
                continue
            for counting_file_row, comparison_file_row in zip(
                csv.reader(open(
                    self._gene_wise_quanti_combined_path), delimiter="\t"),
                    csv.reader(deseq2_result_fh, delimiter="\t")):
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

    def create_final_output_files(self):
        for output_file in os.listdir(self._deseq2_extended_folder):
            comments_csv = pd.read_csv('{}/{}'.format(
                self._deseq2_extended_folder, output_file), sep='\t', nrows=1)
            data_csv = pd.read_csv('{}/{}'.format(
                self._deseq2_extended_folder, output_file), sep='\t', header=2)
            data_csv['FoldChange'] = 2 ** data_csv['log2FoldChange']
            for lib in self._libs:
                if lib in data_csv.columns:
                    data_csv = data_csv.rename(columns={
                        lib: '{} countings'.format(lib)})
            new_output_csv = open('{}/{}'.format(
                self._deseq2_extended_folder, output_file), 'a')
            new_output_csv.write('{}\n'.format(''.join(list(comments_csv))))
            new_output_csv.write('{}\n'.format(''.join(list(
                comments_csv.iloc[0]))))
            data_csv.to_csv(new_output_csv, sep='\t')
            new_output_csv.close()
            # os.remove('{}/{}'.format(self._deseq2_extended_folder, output_file))
        
    def _condition_combos(self, conditions):
        non_redundant_conditions = set(conditions)
        for cond1 in non_redundant_conditions:
            for cond2 in non_redundant_conditions:
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
            comparison_file = "DESeq2_comp_%s_vs_%s.csv" % (
                condition_combo[0], condition_combo[1])
            self._comparison_files_and_combos.append(
                (comparison_file, list(condition_combo)))
            call_string += (
                "write.table(comp%s, file='%s/%s', "
                "quote=FALSE, sep='\\t')\n" % (
                    index, self._deseq2_raw_folder, comparison_file))
        return call_string

    def _deseq2_script_template(self):
        return (
            "library('DESeq2')\n"
            "rawCountTable <- read.table('%s', skip=1, sep='\\t', "
            "quote='', comment.char='', "
            "colClasses=c(rep('character',%s), rep('numeric',%s)))\n"
            "countTable <- round(rawCountTable[,%s:length(names("
            "rawCountTable))])\n"
            "libs <- c(%s)\n"
            "conds <- c(%s)\n"
            "colnames(countTable) <- libs\n"
            "samples <- data.frame(row.names=libs, condition=conds, "
            "lib=libs)\n"
            "dds <- DESeqDataSetFromMatrix(countData=countTable, "
            "colData=samples, design=~condition)\n"
            "dds <- DESeq(dds)\n\n"
            "# PCA plot\n"
            "pdf('%s')\n"
            "rld <- rlog(dds)\n"
            "print(plotPCA(rld, intgroup=c('condition')))\n\n"
            "# Heatmap\n"
            "library('RColorBrewer')\n"
            "library('gplots')\n"
            "distsRL <- dist(t(assay(rld)))\n"
            "mat <- as.matrix(distsRL)\n"
            "rownames(mat) <- with(colData(dds), "
            "paste(lib, sep=' : '))\n"
            "hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)\n"
            "heatmap.2(mat, trace='none', col = rev(hmcol), "
            "margin=c(13, 13))\n")
