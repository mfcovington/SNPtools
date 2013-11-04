
geno_df <- data.frame()
for (i in 1:length(filenames)) {
    # check if (file.info(filenames[i])$size > 0)
    geno_df <- rbind(geno_df,read.table(filenames[i]))
}
colnames(geno_df) <- c("chr","pos","par1", "par2", "cov")

geno_df$cov.plot <- geno_df$cov
geno_df$cov.plot[geno_df$par2 - geno_df$par1 < 0] <- -geno_df$cov.plot[geno_df$par2 - geno_df$par1 < 0]


