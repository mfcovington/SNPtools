# TODO:
# - change chromosome

IL_filenames <- list.files(pattern = "IL2-.*\\.SL2\\.40ch02\\.genotyped\\.nr")
chr          <- "SL2.40ch02"
IL           <- sapply(strsplit(IL_filenames, "\\."), "[", 1)
data_chr     <- as.numeric(substr(chr, 9, 10))
IL_chr       <- sapply(strsplit(IL, "IL|-"), "[", 2)
IL_filenames <- cbind(IL_filenames, IL, data_chr, IL_chr)
IL_filenames <- as.data.frame(IL_filenames)
colnames(IL_filenames) <- c("filename", "IL", "data_chr", "IL_chr")

valid_ILs <- unique(IL_filenames$IL)


### Plot All relevant chromosome, 1 chr (all relevant ILs) per plot
offset <- 1.05
chromosome <- sort(unique(IL_filenames$IL_chr))  # change this


ILs_to_input <- as.character(
  IL_filenames[
    (IL_filenames$IL_chr == chromosome) &
    (IL_filenames$IL_chr == IL_filenames$data_chr)
  , ]$filename
)

IL_df <- data.frame()
for (j in 1:length(ILs_to_input)) {
    if (file.info(ILs_to_input[j])$size > 0) {
        IL_df_temp <- read.table(ILs_to_input[j])
        IL_df_temp <- cbind(IL_df_temp, sapply(strsplit(ILs_to_input[j], "\\."), "[", 1))
        IL_df <- rbind(IL_df,IL_df_temp)
    }
}
colnames(IL_df) <- c("chr","pos","par1", "par2", "cov", "IL")

IL_df$cov.plot <- IL_df$cov
IL_df$cov.plot[(IL_df$par2 - IL_df$par1) / IL_df$cov < 0 & IL_df$cov > 0] <- -IL_df$cov.plot[(IL_df$par2 - IL_df$par1) / IL_df$cov < 0 & IL_df$cov > 0] #CHANGED
chr_formatted <- paste("chr", substr(IL_df$chr,9,10), sep = "")
IL_df$chr <- chr_formatted

