library(DESeq2)

dds <- DESeqDataSetFromMatrix(expreset_agg_round, DataFrame(condition), design= ~ condition )
head(dds)
dds <- DESeq(dds)
res <- results(dds)
summary(res)
gene_symbol <- read.csv("Data\\gene_symbol_all_22268.csv")
sig_id <- read.csv("sigid_A549.csv",header=FALSE)
