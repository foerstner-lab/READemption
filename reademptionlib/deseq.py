import csv
import sys
import os
from subprocess import call


class DESeqRunner(object):
    def __init__(
        self,
        libs,
        conditions,
        deseq_raw_folder,
        deseq_extended_folder,
        deseq_script_path,
        deseq_pca_heatmap_path,
        gene_wise_quanti_combined_path,
        deseq_tmp_session_info_script,
        deseq_session_info,
        fc_shrinkage_off,
        cooks_cutoff_off=False,
    ):
        self._libs = libs
        self._conditions = conditions
        self._deseq_raw_folder = deseq_raw_folder
        self._deseq_extended_folder = deseq_extended_folder
        self._deseq_script_path = deseq_script_path
        self._deseq_pca_heatmap_path = deseq_pca_heatmap_path
        self._gene_wise_quanti_combined_path = gene_wise_quanti_combined_path
        self._deseq_tmp_session_info_script = deseq_tmp_session_info_script
        self._deseq_session_info = deseq_session_info
        self._cooks_cutoff_off = cooks_cutoff_off
        self._first_data_column = 11
        self._fc_shrinkage_off = fc_shrinkage_off

    def write_session_info_file(self):
        with open(self._deseq_tmp_session_info_script, "w") as tmp_r_script_fh:
            tmp_r_script_fh.write("library('DESeq2')\nsessionInfo()\n")
        with open(self._deseq_session_info, "w") as session_info_fh:
            with open(os.devnull, "w") as devnull:
                call(
                    ["Rscript", self._deseq_tmp_session_info_script],
                    stdout=session_info_fh,
                    stderr=devnull,
                )
        os.remove(self._deseq_tmp_session_info_script)

    def create_deseq_script_file(self):
        #{'Co_culture_replicate_1': 'co_culture', 'Co_culture_replicate_2': 'co_culture', 'Co_culture_replicate_3': 'co_culture', 'Harboring_replicate_1': 'harboring', 'Harboring_replicate_2': 'harboring', 'Harboring_replicate_3': 'harboring', 'Infected_replicate_1': 'infected', 'Infected_replicate_2': 'infected', 'Infected_replicate_3': 'infected', 'Steady_state_replicate_1': 'steady_state', 'Steady_state_replicate_2': 'steady_state', 'Steady_state_replicate_3': 'steady_state', 'Uninfected_replicate_1': 'uninfected', 'Uninfected_replicate_2': 'uninfected', 'Uninfected_replicate_3': 'uninfected'}
        libs_to_conditions = dict(
            [
                (lib, condition)
                for lib, condition in zip(self._libs, self._conditions)
            ]
        )
        #['Orientation of counted reads relative to the strand location of the annotation', 'Sequence name', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attributes', 'Co_culture_replicate_1', 'Co_culture_replicate_2', 'Co_culture_replicate_3', 'Harboring_replicate_1', 'Harboring_replicate_2', 'Harboring_replicate_3', 'Infected_replicate_1', 'Infected_replicate_2', 'Infected_replicate_3', 'Steady_state_replicate_1', 'Steady_state_replicate_2', 'Steady_state_replicate_3', 'Uninfected_replicate_1', 'Uninfected_replicate_2', 'Uninfected_replicate_3']
        head_row = (
            open(self._gene_wise_quanti_combined_path)
            .readline()[:-1]
            .split("\t")
        )
        #['Co_culture_replicate_1', 'Co_culture_replicate_2', 'Co_culture_replicate_3', 'Harboring_replicate_1', 'Harboring_replicate_2', 'Harboring_replicate_3', 'Infected_replicate_1', 'Infected_replicate_2', 'Infected_replicate_3', 'Steady_state_replicate_1', 'Steady_state_replicate_2', 'Steady_state_replicate_3', 'Uninfected_replicate_1', 'Uninfected_replicate_2', 'Uninfected_replicate_3']
        libs = head_row[self._first_data_column - 1 :]
        #'Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3','Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3','Infected_replicate_1','Infected_replicate_2','Infected_replicate_3','Steady_state_replicate_1','Steady_state_replicate_2','Steady_state_replicate_3','Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3'
        libs_str = ",".join([f"'{lib}'" for lib in libs])
        conditions = [libs_to_conditions[lib] for lib in libs]
        condition_str = ", ".join([f"'{cond}'" for cond in conditions])
        if not self._fc_shrinkage_off:
            beta_prior_str = ", betaPrior=TRUE"
        elif self._fc_shrinkage_off:
            beta_prior_str = ""
        file_content = self._deseq_script_template(len(libs), libs_str, condition_str, beta_prior_str)
        file_content += self._comparison_call_strings(conditions)
        deseq_fh = open(self._deseq_script_path, "w")
        deseq_fh.write(file_content)
        deseq_fh.close()

    def run_deseq(self):
        call(["Rscript", self._deseq_script_path])

    def merge_counting_files_with_results(self):
        for comparison_file, combo in self._comparison_files_and_combos:

            output_fh = open(f"{self._deseq_extended_folder}/{comparison_file.replace('.csv', '_with_annotation_and_countings.csv')}", "w")
            output_fh.write(
                f"# Reference library (divisor): {combo[1]}\n"
            )
            output_fh.write(
                f"# Comparison library (numerator): {combo[0]}\n"
            )
            try:
                deseq_result_fh = open(
                    f"{self._deseq_raw_folder}/{comparison_file}"
                )
            except:
                sys.stderr.write(
                    f"Apparently DESeq did not generate the "
                    f"file '{comparison_file}'. Extension stopped.\n"
                )
                continue
            for counting_file_row, comparison_file_row in zip(
                csv.reader(
                    open(self._gene_wise_quanti_combined_path), delimiter="\t"
                ),
                csv.reader(deseq_result_fh, delimiter="\t"),
            ):
                if comparison_file_row[0] == "baseMean":
                    # Add another column to the header
                    comparison_file_row = [""] + comparison_file_row
                    # Extend column description
                    counting_file_row[self._first_data_column :] = [
                        f"{lib_name} raw countings"
                        for lib_name in counting_file_row[
                            self._first_data_column :
                        ]
                    ]
                output_fh.write(
                    "\t".join(counting_file_row + comparison_file_row[1:])
                    + "\n"
                )
            output_fh.close()

    def _condition_combos(self, conditions):
        non_redundant_conditions = set(conditions)
        for cond1 in non_redundant_conditions:
            for cond2 in non_redundant_conditions:
                if not cond1 == cond2:
                    yield ((cond1, cond2))

    def _comparison_call_strings(self, conditions):
        call_string = ""
        condition_combos = self._condition_combos(conditions)
        self._comparison_files_and_combos = []
        cooks_cutoff_str = ""
        if self._cooks_cutoff_off:
            cooks_cutoff_str = ", cooksCutoff=FALSE"
        for index, condition_combo in enumerate(condition_combos):
            call_string += (
                f"comp{index} <- results(dds, contrast="
                f"c('condition','{condition_combo[0]}', '{condition_combo[1]}'){cooks_cutoff_str})\n"
            )
            comparison_file = f"deseq_comp_{condition_combo[0]}_vs_{condition_combo[1]}.csv"
            self._comparison_files_and_combos.append(
                (comparison_file, list(condition_combo))
            )
            call_string += (
                f"write.table(comp{index}, file='{self._deseq_raw_folder}/{comparison_file}', "
                "quote=FALSE, sep='\\t')\n"
            )
        return call_string


    #self._gene_wise_quanti_combined_path,
    #self._first_data_column - 1,
    #len(libs),
    #self._first_data_column,
    #libs_str,
    #condition_str,
    #beta_prior_str,
    #self._deseq_pca_heatmap_path

    def _deseq_script_template(self, libs_length, libs_str, condition_str, beta_prior_str):
        return (
            f"library('DESeq2')\n"
            f"rawCountTable <- read.table('{self._gene_wise_quanti_combined_path}', skip=1, sep='\\t', "
            f"quote='', comment.char='', "
            f"colClasses=c(rep('character',{self._first_data_column - 1}), rep('numeric',{libs_length})))\n"
            f"countTable <- round(rawCountTable[,{self._first_data_column}:length(names("
            f"rawCountTable))])\n"
            f"libs <- c({libs_str})\n"
            f"conds <- c({condition_str})\n"
            f"colnames(countTable) <- libs\n"
            f"samples <- data.frame(row.names=libs, condition=conds, "
            f"lib=libs)\n"
            f"head(samples)\n"
            f"dds <- DESeqDataSetFromMatrix(countData=countTable, "
            f"colData=samples, design=~condition)\n"
            f"dds <- DESeq(dds{beta_prior_str})\n\n"
            f"# PCA plot\n"
            f"pdf('{self._deseq_pca_heatmap_path}')\n"
            f"rld <- rlog(dds)\n"
            f"print(plotPCA(rld, intgroup=c('condition')))\n\n"
            f"# Heatmap\n"
            f"library('RColorBrewer')\n"
            f"library('gplots')\n"
            f"distsRL <- dist(t(assay(rld)))\n"
            f"mat <- as.matrix(distsRL)\n"
            f"rownames(mat) <- with(colData(dds), "
            f"paste(lib, sep=' : '))\n"
            f"hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)\n"
            f"heatmap.2(mat, trace='none', col = rev(hmcol), "
            f"margin=c(13, 13))\n"
        )
