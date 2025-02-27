#!/usr/bin/env Rscript
# File     :  calSIG.R
# Time     :  2024/10/01 08:08:00
# Author   :  Wenyong Zhu
# Version  :  1.2.0
# Desc     :  calculate the consistence of transcriptional signal across the sample(s)

rm(list = ls())

args <- commandArgs(TRUE)
sample_label <- args[1]
list_bigwig <- args[2]

bw_files <- read.table(list_bigwig, header = FALSE, col.names = c("bw"), stringsAsFactors = FALSE)
total_num <- nrow(bw_files)

t_sig <- data.frame()
p_val <- data.frame()
for (bw_file in bw_files$bw){
    # random transcriptional signal
    df <- read.table(paste0(bw_file, ".rdsig"), header = FALSE)[, 2]
    # average transcriptional signal from bigWigAverageOverBed
    fn <- ecdf(df)
    signal <- read.table(paste0(bw_file, ".avgsig"), header = FALSE)
    pvalue <- as.numeric(format(1-fn(signal[,2]), digits = 3))
    # merge
    if (ncol(t_sig) == 0) {t_sig <- signal; signal[,2] <- pvalue; p_val <- signal;}
    else {t_sig <- cbind(t_sig, signal[,2]); p_val <- cbind(p_val, pvalue);}
}
colnames(t_sig) <- c("locus", sub(".*/([^/]+)\\.bw", "\\1", bw_files$bw)); colnames(p_val) <- c("locus", sub(".*/([^/]+)\\.bw", "\\1", bw_files$bw))
write.table(t_sig, paste0("TNE.", sample_label, ".avgsig.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(p_val, paste0("TNE.", sample_label, ".pvalues.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

## binomial test for the significant TNE (p<=0.05)
if (total_num == 1) {
    binomial_pvalue <- sapply(as.numeric(p_val[, -1] <= 0.05), function(x) binom.test(x, 1, 0.05, 'greater')$p.value)
} else {
    binomial_pvalue <- sapply(rowSums(p_val[, -1] <= 0.05), function(x) binom.test(x, total_num, 0.05, 'greater')$p.value)
}

# using the Holm-Bonferroni method (AKA step-down Bonferroni) to correct for multiple test
adjusted_p <- cbind(binomial.pvalues = binomial_pvalue, p.adjusted.HB = p.adjust(binomial_pvalue, method = "holm"), p.adjusted.bonferroni = p.adjust(binomial_pvalue, method = "bonferroni"), p.adjusted.FDR = p.adjust(binomial_pvalue, method = "fdr"))
rownames(adjusted_p) <- p_val[, 1]
write.table(adjusted_p, paste0("TNE.", sample_label, ".adj_pvalues.txt"), col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)