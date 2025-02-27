#!/usr/bin/env Rscript
# File     :  fitRTS.R
# Time     :  2024/10/01 08:08:00
# Author   :  Wenyong Zhu
# Version  :  1.2.0
# Desc     :  fit distribution of random trancriptional signal

rm(list = ls())

suppressPackageStartupMessages(library(fitdistrplus))

args <- commandArgs(TRUE)
transcriptional_noise <- args[1]

df <- read.table(transcriptional_noise, comment.char = "")
df <- log10(df[, 1] + 1e-16)

fitn <- fitdist(df, "norm")

m <- round(as.numeric(fitn$estimate[1]), digits = 3)
sd <- round(as.numeric(fitn$estimate[2]), digits = 3)

p <- round(qnorm(0.05, mean = m, sd = sd, lower.tail = FALSE), digits = 3)

write.table(10**p - 1e-16, "transcriptional_noise_pvalues.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)