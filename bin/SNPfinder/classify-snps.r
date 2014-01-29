#!/usr/bin/env Rscript --vanilla

args <- commandArgs(trailingOnly = TRUE)
work.dir <- args[1]
setwd(work.dir)
library("foreach")

classify_snps <- function(filename_pattern) {
  allDup <- function (value) { duplicated(value) | duplicated(value, fromLast = TRUE) }

  foreach(file_name = list.files(pattern = filename_pattern)) %do% {
    snps <- read.table(file_name, header = T)
    snps <- cbind(snps,"SNP_CLASS" = factor(NA, levels = c("SNP", "DIFF_SNP", "NOT", "CONFLICT")))
    try({snps[!allDup(snps[,c(1:2,6)]),]$SNP_CLASS = "SNP"}, silent = T)
    try({snps[as.logical(allDup(snps[,c(1,2,6)]) - allDup(snps[,c(1:4,6)])), ]$SNP_CLASS = "DIFF_SNP"}, silent = T)
    try({snps[as.logical(allDup(snps[,c(1:4,6)]) - allDup(snps[,1:6])), ]$SNP_CLASS = "NOT"}, silent = T)
    try({snps[as.logical(allDup(snps[,c(1,2,5:6)]) - allDup(snps[,1:6])), ]$SNP_CLASS = "CONFLICT"}, silent = T)
    print(summary(snps))
    write.table(snps, file = paste(file_name, ".classified", sep = ""), quote = F, row.names = F, sep = "\t")
  }
}

classify_snps("master_snp_list")
