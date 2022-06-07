library('DESeq2')
library('RColorBrewer')
library('gplots')
library('ggplot2')
rawCountTable <- read.table('a_test_project/output/human_gene_quanti_combined/gene_wise_quantifications_combined.csv', skip=1, sep='\t', quote='', comment.char='', colClasses=c(rep('character',10), rep('numeric',15)))
countTable <- round(rawCountTable[,11:length(names(rawCountTable))])
colnames(countTable) <- c('Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3','Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3','Infected_replicate_1','Infected_replicate_2','Infected_replicate_3','Steady_state_replicate_1','Steady_state_replicate_2','Steady_state_replicate_3','Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3')
# Select only the libraries of this species
countTable <- countTable[, c('Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3','Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3','Infected_replicate_1','Infected_replicate_2','Infected_replicate_3','Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3')]
libs <- c('Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3','Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3','Infected_replicate_1','Infected_replicate_2','Infected_replicate_3','Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3')
conds <- c('co_culture', 'co_culture', 'co_culture', 'harboring', 'harboring', 'harboring', 'infected', 'infected', 'infected', 'uninfected', 'uninfected', 'uninfected')
reps <- c('1', '2', '3', '1', '2', '3', '1', '2', '3', '1', '2', '3')
samples <- data.frame(row.names=libs, condition=conds, lib=libs, replicate=reps)
dds <- DESeqDataSetFromMatrix(countData=countTable, colData=samples, design=~condition)
dds <- DESeq(dds, betaPrior=TRUE)

# PCA plot
pdf('a_test_project/output/human_deseq/deseq_raw/sample_comparison_pca_heatmap.pdf')
rld <- rlog(dds)
pcaData <- plotPCA(rld, 'condition', intgroup=c('condition', 'replicate'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, 'percentVar'))
print(ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
geom_point(size=3) +
xlab(paste0('PC1: ',percentVar[1],'% variance')) +
ylab(paste0('PC2: ',percentVar[2],'% variance')) +
coord_fixed())
# Heatmap
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- with(colData(dds), paste(lib, sep=' : '))
hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmap.2(mat, trace='none', col = rev(hmcol), margin=c(13, 13))
#Comparison: co_culture vs. harboring 
countTable0 <- countTable[, c('Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3','Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3')]
libs0 <- c('Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3','Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3')
conds0 <- c('co_culture','co_culture','co_culture','harboring','harboring','harboring')
samples0 <- data.frame(row.names=libs0, condition=conds0, lib=libs0)
dds0 <- DESeqDataSetFromMatrix(countData=countTable0, colData=samples0, design=~condition)
dds0 <- DESeq(dds0, betaPrior=TRUE)

comp0 <- results(dds0, contrast=c('condition','co_culture', 'harboring'))
write.table(comp0, file='a_test_project/output/human_deseq/deseq_raw/deseq_comp_co_culture_vs_harboring.csv', quote=FALSE, sep='\t')
#Comparison: co_culture vs. infected 
countTable1 <- countTable[, c('Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3','Infected_replicate_1','Infected_replicate_2','Infected_replicate_3')]
libs1 <- c('Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3','Infected_replicate_1','Infected_replicate_2','Infected_replicate_3')
conds1 <- c('co_culture','co_culture','co_culture','infected','infected','infected')
samples1 <- data.frame(row.names=libs1, condition=conds1, lib=libs1)
dds1 <- DESeqDataSetFromMatrix(countData=countTable1, colData=samples1, design=~condition)
dds1 <- DESeq(dds1, betaPrior=TRUE)

comp1 <- results(dds1, contrast=c('condition','co_culture', 'infected'))
write.table(comp1, file='a_test_project/output/human_deseq/deseq_raw/deseq_comp_co_culture_vs_infected.csv', quote=FALSE, sep='\t')
#Comparison: co_culture vs. uninfected 
countTable2 <- countTable[, c('Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3','Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3')]
libs2 <- c('Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3','Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3')
conds2 <- c('co_culture','co_culture','co_culture','uninfected','uninfected','uninfected')
samples2 <- data.frame(row.names=libs2, condition=conds2, lib=libs2)
dds2 <- DESeqDataSetFromMatrix(countData=countTable2, colData=samples2, design=~condition)
dds2 <- DESeq(dds2, betaPrior=TRUE)

comp2 <- results(dds2, contrast=c('condition','co_culture', 'uninfected'))
write.table(comp2, file='a_test_project/output/human_deseq/deseq_raw/deseq_comp_co_culture_vs_uninfected.csv', quote=FALSE, sep='\t')
#Comparison: harboring vs. co_culture 
countTable3 <- countTable[, c('Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3','Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3')]
libs3 <- c('Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3','Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3')
conds3 <- c('harboring','harboring','harboring','co_culture','co_culture','co_culture')
samples3 <- data.frame(row.names=libs3, condition=conds3, lib=libs3)
dds3 <- DESeqDataSetFromMatrix(countData=countTable3, colData=samples3, design=~condition)
dds3 <- DESeq(dds3, betaPrior=TRUE)

comp3 <- results(dds3, contrast=c('condition','harboring', 'co_culture'))
write.table(comp3, file='a_test_project/output/human_deseq/deseq_raw/deseq_comp_harboring_vs_co_culture.csv', quote=FALSE, sep='\t')
#Comparison: harboring vs. infected 
countTable4 <- countTable[, c('Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3','Infected_replicate_1','Infected_replicate_2','Infected_replicate_3')]
libs4 <- c('Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3','Infected_replicate_1','Infected_replicate_2','Infected_replicate_3')
conds4 <- c('harboring','harboring','harboring','infected','infected','infected')
samples4 <- data.frame(row.names=libs4, condition=conds4, lib=libs4)
dds4 <- DESeqDataSetFromMatrix(countData=countTable4, colData=samples4, design=~condition)
dds4 <- DESeq(dds4, betaPrior=TRUE)

comp4 <- results(dds4, contrast=c('condition','harboring', 'infected'))
write.table(comp4, file='a_test_project/output/human_deseq/deseq_raw/deseq_comp_harboring_vs_infected.csv', quote=FALSE, sep='\t')
#Comparison: harboring vs. uninfected 
countTable5 <- countTable[, c('Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3','Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3')]
libs5 <- c('Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3','Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3')
conds5 <- c('harboring','harboring','harboring','uninfected','uninfected','uninfected')
samples5 <- data.frame(row.names=libs5, condition=conds5, lib=libs5)
dds5 <- DESeqDataSetFromMatrix(countData=countTable5, colData=samples5, design=~condition)
dds5 <- DESeq(dds5, betaPrior=TRUE)

comp5 <- results(dds5, contrast=c('condition','harboring', 'uninfected'))
write.table(comp5, file='a_test_project/output/human_deseq/deseq_raw/deseq_comp_harboring_vs_uninfected.csv', quote=FALSE, sep='\t')
#Comparison: infected vs. co_culture 
countTable6 <- countTable[, c('Infected_replicate_1','Infected_replicate_2','Infected_replicate_3','Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3')]
libs6 <- c('Infected_replicate_1','Infected_replicate_2','Infected_replicate_3','Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3')
conds6 <- c('infected','infected','infected','co_culture','co_culture','co_culture')
samples6 <- data.frame(row.names=libs6, condition=conds6, lib=libs6)
dds6 <- DESeqDataSetFromMatrix(countData=countTable6, colData=samples6, design=~condition)
dds6 <- DESeq(dds6, betaPrior=TRUE)

comp6 <- results(dds6, contrast=c('condition','infected', 'co_culture'))
write.table(comp6, file='a_test_project/output/human_deseq/deseq_raw/deseq_comp_infected_vs_co_culture.csv', quote=FALSE, sep='\t')
#Comparison: infected vs. harboring 
countTable7 <- countTable[, c('Infected_replicate_1','Infected_replicate_2','Infected_replicate_3','Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3')]
libs7 <- c('Infected_replicate_1','Infected_replicate_2','Infected_replicate_3','Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3')
conds7 <- c('infected','infected','infected','harboring','harboring','harboring')
samples7 <- data.frame(row.names=libs7, condition=conds7, lib=libs7)
dds7 <- DESeqDataSetFromMatrix(countData=countTable7, colData=samples7, design=~condition)
dds7 <- DESeq(dds7, betaPrior=TRUE)

comp7 <- results(dds7, contrast=c('condition','infected', 'harboring'))
write.table(comp7, file='a_test_project/output/human_deseq/deseq_raw/deseq_comp_infected_vs_harboring.csv', quote=FALSE, sep='\t')
#Comparison: infected vs. uninfected 
countTable8 <- countTable[, c('Infected_replicate_1','Infected_replicate_2','Infected_replicate_3','Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3')]
libs8 <- c('Infected_replicate_1','Infected_replicate_2','Infected_replicate_3','Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3')
conds8 <- c('infected','infected','infected','uninfected','uninfected','uninfected')
samples8 <- data.frame(row.names=libs8, condition=conds8, lib=libs8)
dds8 <- DESeqDataSetFromMatrix(countData=countTable8, colData=samples8, design=~condition)
dds8 <- DESeq(dds8, betaPrior=TRUE)

comp8 <- results(dds8, contrast=c('condition','infected', 'uninfected'))
write.table(comp8, file='a_test_project/output/human_deseq/deseq_raw/deseq_comp_infected_vs_uninfected.csv', quote=FALSE, sep='\t')
#Comparison: uninfected vs. co_culture 
countTable9 <- countTable[, c('Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3','Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3')]
libs9 <- c('Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3','Co_culture_replicate_1','Co_culture_replicate_2','Co_culture_replicate_3')
conds9 <- c('uninfected','uninfected','uninfected','co_culture','co_culture','co_culture')
samples9 <- data.frame(row.names=libs9, condition=conds9, lib=libs9)
dds9 <- DESeqDataSetFromMatrix(countData=countTable9, colData=samples9, design=~condition)
dds9 <- DESeq(dds9, betaPrior=TRUE)

comp9 <- results(dds9, contrast=c('condition','uninfected', 'co_culture'))
write.table(comp9, file='a_test_project/output/human_deseq/deseq_raw/deseq_comp_uninfected_vs_co_culture.csv', quote=FALSE, sep='\t')
#Comparison: uninfected vs. harboring 
countTable10 <- countTable[, c('Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3','Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3')]
libs10 <- c('Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3','Harboring_replicate_1','Harboring_replicate_2','Harboring_replicate_3')
conds10 <- c('uninfected','uninfected','uninfected','harboring','harboring','harboring')
samples10 <- data.frame(row.names=libs10, condition=conds10, lib=libs10)
dds10 <- DESeqDataSetFromMatrix(countData=countTable10, colData=samples10, design=~condition)
dds10 <- DESeq(dds10, betaPrior=TRUE)

comp10 <- results(dds10, contrast=c('condition','uninfected', 'harboring'))
write.table(comp10, file='a_test_project/output/human_deseq/deseq_raw/deseq_comp_uninfected_vs_harboring.csv', quote=FALSE, sep='\t')
#Comparison: uninfected vs. infected 
countTable11 <- countTable[, c('Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3','Infected_replicate_1','Infected_replicate_2','Infected_replicate_3')]
libs11 <- c('Uninfected_replicate_1','Uninfected_replicate_2','Uninfected_replicate_3','Infected_replicate_1','Infected_replicate_2','Infected_replicate_3')
conds11 <- c('uninfected','uninfected','uninfected','infected','infected','infected')
samples11 <- data.frame(row.names=libs11, condition=conds11, lib=libs11)
dds11 <- DESeqDataSetFromMatrix(countData=countTable11, colData=samples11, design=~condition)
dds11 <- DESeq(dds11, betaPrior=TRUE)

comp11 <- results(dds11, contrast=c('condition','uninfected', 'infected'))
write.table(comp11, file='a_test_project/output/human_deseq/deseq_raw/deseq_comp_uninfected_vs_infected.csv', quote=FALSE, sep='\t')
