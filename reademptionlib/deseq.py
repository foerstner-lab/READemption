import csv
import sys
import os
from subprocess import call


class DESeqRunner(object):
    def __init__(
        self,
        species,
        libs,
        conditions,
        replicates,
        libs_by_species,
        size_factor,
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
        self._species = species
        self._libs = libs
        self._conditions = conditions
        self._replicates = replicates
        self._libs_by_species = libs_by_species
        self._size_factor = size_factor
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
        libs_to_conditions = dict(
            [
                (lib, condition)
                for lib, condition in zip(self._libs, self._conditions)
            ]
        )
        libs_to_replicates = dict(
            [(lib, rep) for lib, rep in zip(self._libs, self._replicates)]
        )
        libs_by_conditions = {}
        for unique_condition in set(self._conditions):
            libs_by_conditions[unique_condition] = []
        for lib, condition in libs_to_conditions.items():
            libs_by_conditions[condition].append(lib)
        head_row = (
            open(self._gene_wise_quanti_combined_path)
            .readline()[:-1]
            .split("\t")
        )
        libs = head_row[self._first_data_column - 1 :]
        libs_str_quanti_table = ",".join([f"'{lib}'" for lib in libs])
        # Create libs-string, conditions-string and replicates-string in the
        # same order as libs appear in the gene quantification table of the
        # ALL species
        if self._size_factor == "project":
            libs_str = ",".join([f"'{lib}'" for lib in libs])
            conditions = [libs_to_conditions[lib] for lib in libs]
            condition_str = ", ".join([f"'{cond}'" for cond in conditions])
            replicates = [libs_to_replicates[lib] for lib in libs]
            replicates_str = ", ".join(
                [f"'{replicate}'" for replicate in replicates]
            )
        # Create libs-string, conditions-string and replicates-string in the
        # same order as libs appear in the gene quantification table of the
        # ONE species
        if not self._size_factor == "project":
            libs_of_species = []
            conditions = []
            replicates = []
            for lib in libs:
                if lib in self._libs_by_species[self._species]:
                    libs_of_species.append(lib)
                    conditions.append(libs_to_conditions[lib])
                    replicates.append(libs_to_replicates[lib])
            libs_str = ",".join([f"'{lib}'" for lib in libs_of_species])
            condition_str = ", ".join([f"'{cond}'" for cond in conditions])
            replicates_str = ", ".join(
                [f"'{replicate}'" for replicate in replicates]
            )
        if not self._fc_shrinkage_off:
            beta_prior_str = ", betaPrior=TRUE"
        elif self._fc_shrinkage_off:
            beta_prior_str = ""
        file_content = self._deseq_script_template(
            libs_str_quanti_table,
            len(libs),
            libs_str,
            condition_str,
            replicates_str,
            beta_prior_str,
        )
        file_content += self._comparison_call_strings(
            conditions,
            libs_by_conditions,
            libs_to_conditions,
            libs_to_replicates,
            beta_prior_str,
        )
        deseq_fh = open(self._deseq_script_path, "w")
        deseq_fh.write(file_content)
        deseq_fh.close()

    def run_deseq(self):
        call(["Rscript", self._deseq_script_path])

    def merge_counting_files_with_results(self):
        for comparison_file, combo in self._comparison_files_and_combos:
            output_fh = open(
                f"{self._deseq_extended_folder}/{comparison_file.replace('.csv', '_with_annotation_and_countings.csv')}",
                "w",
            )
            output_fh.write(f"# Reference library (divisor): {combo[1]}\n")
            output_fh.write(f"# Comparison library (numerator): {combo[0]}\n")
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
                    counting_file_row[(self._first_data_column -1) :] = [
                        f"{lib_name} raw countings"
                        for lib_name in counting_file_row[
                            (self._first_data_column -1) :
                        ]
                    ]
                output_fh.write(
                    "\t".join(counting_file_row + comparison_file_row[1:])
                    + "\n"
                )
            output_fh.close()

    def _condition_combos(self, conditions):
        non_redundant_conditions = []
        for cond in conditions:
            if cond not in non_redundant_conditions:
                non_redundant_conditions.append(cond)
        for cond1 in non_redundant_conditions:
            for cond2 in non_redundant_conditions:
                if not cond1 == cond2:
                    yield ((cond1, cond2))

    def _comparison_call_strings(
        self,
        conditions,
        libs_by_conditions,
        libs_to_conditions,
        libs_to_replicates,
        beta_prior_str,
    ):
        call_string = ""
        condition_combos = self._condition_combos(conditions)
        self._comparison_files_and_combos = []
        cooks_cutoff_str = ""
        if self._cooks_cutoff_off:
            cooks_cutoff_str = ", cooksCutoff=FALSE"
        for index, condition_combo in enumerate(condition_combos):
            if self._size_factor == "comparison":
                # Add the single comparisons based only on the libs being
                # compared (a dds for each comparison)
                lib_names_to_keep = []
                for condition in condition_combo:
                    lib_names_to_keep.extend(libs_by_conditions[condition])
                lib_names_to_keep_string = ",".join(
                    [f"'{lib_name}'" for lib_name in lib_names_to_keep]
                )
                conditions_to_keep = []
                for lib in lib_names_to_keep:
                    conditions_to_keep.append(libs_to_conditions[lib])
                conditions_to_keep_string = ",".join(
                    [f"'{condition}'" for condition in conditions_to_keep]
                )
                replicates_to_keep = []
                for lib in lib_names_to_keep:
                    replicates_to_keep.append(libs_to_replicates[lib])
                replicates_to_keep_string = ",".join(
                    [f"'{replicate}'" for replicate in replicates_to_keep]
                )
                call_string += (
                    f"#Comparison: {condition_combo[0]} vs. {condition_combo[1]} \n"
                    # Select only the libraries to compare and create a new dataframe
                    f"countTable{index} <- countTable[, c({lib_names_to_keep_string})]\n"
                    f"libs{index} <- c({lib_names_to_keep_string})\n"
                    f"conds{index} <- c({conditions_to_keep_string})\n"
                    f"samples{index} <- data.frame(row.names=libs{index}, condition=conds{index}, "
                    f"lib=libs{index})\n"
                    f"dds{index} <- DESeqDataSetFromMatrix(countData=countTable{index}, "
                    f"colData=samples{index}, design=~condition)\n"
                    f"dds{index} <- DESeq(dds{index}{beta_prior_str})\n\n"
                )
                # Calculate the DESeq-Data-Set only for the libs of the current
                # comparison (but only for the current species)
                call_string += (
                    f"comp{index} <- results(dds{index}, contrast="
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

            if self._size_factor in ["project", "species"]:
                # Add the single comparisons based on the whole project (dds)
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

    def _deseq_script_template(
        self,
        libs_str_quanti_table,
        libs_length,
        libs_str,
        condition_str,
        replicates_str,
        beta_prior_str,
    ):
        call_string = (
            # Load packages
            f"library('DESeq2')\n"
            f"library('RColorBrewer')\n"
            f"library('gplots')\n"
            f"library('ggplot2')\n"
            # Read in the gene quanti table for the current species
            f"rawCountTable <- read.table('{self._gene_wise_quanti_combined_path}', skip=1, sep='\\t', "
            f"quote='', comment.char='', "
            f"colClasses=c(rep('character',{self._first_data_column - 1}), rep('numeric',{libs_length})))\n"
            f"countTable <- round(rawCountTable[,{self._first_data_column}:length(names("
            f"rawCountTable))])\n"
            f"colnames(countTable) <- c({libs_str_quanti_table})\n"
        )
        if not self._size_factor == "project":
            call_string += (
                f"# Select only the libraries of this species\n"
                f"countTable <- countTable[, c({libs_str})]\n"
            )
        call_string += (
            # Set the libs, conds and reps of the whole project
            f"libs <- c({libs_str})\n"
            f"conds <- c({condition_str})\n"
            f"reps <- c({replicates_str})\n"
            f"samples <- data.frame(row.names=libs, condition=conds, "
            f"lib=libs, replicate=reps)\n"
            # Calculate the DESeq-Data-Set for all libs of the project (but only
            # for the current species)
            f"dds <- DESeqDataSetFromMatrix(countData=countTable, "
            f"colData=samples, design=~condition)\n"
            f"dds <- DESeq(dds{beta_prior_str})\n\n"
        )

        call_string += (
            # Create a PCA plot for all libs of the project (but only
            # for the current species)
            f"# PCA plot\n"
            f"pdf('{self._deseq_pca_heatmap_path}')\n"
            f"rld <- rlog(dds)\n"
            f"pcaData <- plotPCA(rld, 'condition', intgroup=c('condition', 'replicate'), returnData=TRUE)\n"
            f"percentVar <- round(100 * attr(pcaData, 'percentVar'))\n"
            f"print(ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +\n"
            f"geom_point(size=3) +\n"
            f"xlab(paste0('PC1: ',percentVar[1],'% variance')) +\n"
            f"ylab(paste0('PC2: ',percentVar[2],'% variance')) +\n"
            f"coord_fixed())\n"
            # Create a Heatmap plot for all libs of the project (but only
            # for the current species)
            f"# Heatmap\n"
            f"distsRL <- dist(t(assay(rld)))\n"
            f"mat <- as.matrix(distsRL)\n"
            f"rownames(mat) <- with(colData(dds), "
            f"paste(lib, sep=' : '))\n"
            f"hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)\n"
            f"heatmap.2(mat, trace='none', col = rev(hmcol), "
            f"margin=c(13, 13))\n"
        )
        return call_string
