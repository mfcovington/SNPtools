##2012-08-26: for plotting final versions of RNAseq

library("ggplot2")

# TODO:
# - add color customization
# - add plot size customization
# - add plot format customization
# - change chromosome
# - 

IL_filenames <- list.files(pattern = "IL2-.*\\.SL2\\.40ch02\\.genotyped\\.nr")
chr          <- "SL2.40ch02"
IL           <- sapply(strsplit(IL_filenames, "\\."), "[", 1)
data_chr     <- as.numeric(substr(chr, 9, 10))
IL_chr       <- sapply(strsplit(IL, "IL|-"), "[", 2)
IL_filenames <- cbind(IL_filenames, IL, data_chr, IL_chr)
IL_filenames <- as.data.frame(IL_filenames)
colnames(IL_filenames) <- c("filename", "IL", "data_chr", "IL_chr")

valid_ILs <- unique(IL_filenames$IL)

boundaries <- read.table("IL_boundaries", head = T)

# > head (boundaries)
#   chr      IL    start      end PEN
# 1   1   IL1-1        1   303076  -1
# 2   1   IL1-1   309021 77720421   1
# 3   1   IL1-1 77721171 90304244  -1
# 4   1 IL1-1-2        1   303076  -1
# 5   1 IL1-1-2   309061  1137671   1
# 6   1 IL1-1-2  1151331 90304244  -1

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



max_cov <- max(IL_df$cov)
ylimits <- max_cov * offset

# added due to broken feature in ggplot 0.9.1:
max_pos  <- max(IL_df$pos)
temp_max <- floor(max_pos / 10000000)
max_pos_fixed <- temp_max * 10000000
label_max <- temp_max * 10

geno.plot <- ggplot(IL_df, aes(pos, cov.plot)) +
  geom_point(size = 1, aes(color = (par2 - par1) / cov)) + 
  facet_grid(IL ~ .) + 
  scale_color_gradient2(
    low = "magenta",
    mid = "black",
    high = "green",
    name = expression((par2 - par1) / total)
  ) +
  guides(color = guide_legend(reverse = TRUE)) + 
  opts(
    title = paste(sep = " ", "Chromosome", chromosome),
    panel.grid.minor = theme_blank(),
    panel.grid.major = theme_blank()
  ) + 
  scale_x_continuous(
    'Position on chromosome (Mb)',
    breaks = seq(0, max_pos_fixed, 10000000),
    labels = seq(0, label_max, 10)
  ) + 
  scale_y_continuous(
    'Coverage',
    breaks = c(-100, 0, 100),
    labels = c("", "", ""),
    limits = c(-ylimits, ylimits)
  ) + 
  geom_segment(
    data   = boundaries[boundaries$chr == chromosomes[i] & boundaries$IL %in% valid_ILs, ],
    facets = IL ~ .,
    color  = "gray50",
    aes(
      x    = start,
      xend = end,
      y    = PEN * ylimits,
      yend = PEN * ylimits)
  )

geno.plot
# ggsave(paste("example.", chr_formatted, ".png", sep = ""), plot = geno.plot, width = 6, height = 5) #CHANGED


