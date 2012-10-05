##2012-08-26: for plotting final versions of RNAseq

library("ggplot2")

#building data frame (of IL filenames, IL IDs, chr of file, chr of IL) and list of IL IDs
# IL_filenames <- system("ls P*.genotyped.chr*", intern = TRUE)
IL_filenames <- system("ls *.genotyped.chr*.nr", intern = TRUE)

#RNAseq:
IL_filenames <- cbind(IL_filenames, matrix(unlist(strsplit(IL_filenames, ".genotyped.")), ncol = 2, byrow = T), unlist(strsplit(IL_filenames, "\\..*"))[unlist(strsplit(IL_filenames, "\\..*")) != ""])
#RESCAN:
#IL_filenames <- cbind(IL_filenames, matrix(unlist(strsplit(IL_filenames, ".genotyped.")), ncol = 2, byrow = T), substring(unlist(strsplit(IL_filenames, "-.*"))[unlist(strsplit(IL_filenames, "-.*")) != ""], 3))

colnames(IL_filenames) <- c("filename", "IL", "data_chr", "IL_chr")
IL_filenames[,4] <- as.numeric(IL_filenames[,4])
#IL_filenames[,4] <- as.numeric(substr(IL_filenames[,4]))
IL_filenames <- as.data.frame(IL_filenames)

IL_filenames$data_chr <- substr(IL_filenames$data_chr,1,5)


#this is to add IL9-3-1 onto the plot for chr 12
IL_filenames$IL_chr[936] = 12


valid_ILs <- unique(paste(sep="", "IL", gsub("\\.", "-", IL_filenames$IL)))


# # IL_filenames$IL_chr <- as.numeric(IL_filenames$IL_chr)
# IL_list <- unique(IL_filenames$IL)

# ### Plot All, 1 IL (12 chr) per plot
# for (i in 1:length(IL_list)) {
#     ILs_to_input <- as.character(IL_filenames[(IL_filenames$IL == IL_list[i]) & (IL_filenames$data_chr != "chr00"), ]$filename)
# #     IL_id <- as.character(IL_filenames$IL[IL_filenames$IL == IL_list[i]][1])
    
#     IL_df <- data.frame()
#     for (j in 1:length(ILs_to_input)) {
#         if (file.info(ILs_to_input[j])$size > 0) IL_df <- rbind(IL_df,read.table(ILs_to_input[j]))
#     }
    
#     if (nrow(IL_df) > 0) {
#         colnames(IL_df) <- c("chr","pos","M", "P", "COV")
#         IL_df$COV.plot <- IL_df$COV #CHANGED
#         IL_df$COV.plot[(IL_df$P-IL_df$M)/IL_df$COV < 0 & IL_df$COV > 0] <- -IL_df$COV.plot[(IL_df$P-IL_df$M)/IL_df$COV < 0 & IL_df$COV > 0] #CHANGED
#         IL_df$chr<-paste("chr", substr(IL_df$chr,9,10), sep = "")
#         max_COV <- max(IL_df$COV) #CHANGED
#         geno <- ggplot(IL_df, aes(pos, COV.plot))
#         geno <- geno + geom_point(size = 1, aes(color = (P-M)/COV)) + scale_color_gradient2(low = "magenta", mid = "black", high = "green", name = expression((Pen - M82)/total)) + facet_grid(chr ~ .)
#         geno <- geno + opts(title = IL_list[i])
#         geno <- geno + scale_x_continuous('Position on chromosome (Mb)', breaks = seq(0, 80000000, 10000000), labels = seq(0, 80, 10)) + scale_y_continuous('Coverage', breaks = c(-100, 0, 100), labels = c("", "", ""), limits = c(-max_COV, max_COV)) #CHANGED
#         geno <- geno + opts(panel.grid.minor=theme_blank(), panel.grid.major=theme_blank())
#         ggsave(paste("../plot_nr/split.", IL_list[i], ".all.pdf", sep = ""), plot = geno, width = 10.5, height = 8) #CHANGED
#         ggsave(paste("../plot_nr/split.", IL_list[i], ".all.png", sep = ""), plot = geno, width = 10.5, height = 8) #CHANGED
#         # ggsave(paste("../plot/split.", IL_list[i], ".all.pdf", sep = ""), plot = geno, width = 10.5, height = 8) #CHANGED
#         # ggsave(paste("../plot/split.", IL_list[i], ".all.png", sep = ""), plot = geno, width = 10.5, height = 8) #CHANGED
#     }
# }


boundaries <- read.table("IL_boundaries", head = T)

# > head (boundaries)
#   chr      IL    start      end PEN
# 1   1   IL1-1        1   303076  -1
# 2   1   IL1-1   309021 77720421   1
# 3   1   IL1-1 77721171 90304244  -1
# 4   1 IL1-1-2        1   303076  -1
# 5   1 IL1-1-2   309061  1137671   1
# 6   1 IL1-1-2  1151331 90304244  -1

### Plot All relevant chromosomes, 1 chr (all relevant ILs) per plot
offset <- 1.05
chromosomes <- sort(unique(IL_filenames$IL_chr))
for (i in 1:length(chromosomes)) { #length(chromosomes)

    # boundaries_chr <- boundaries[boundaries$chr == chromosomes,]

    ILs_to_input <- as.character(IL_filenames[(IL_filenames$IL_chr == chromosomes[i]) & (IL_filenames$IL_chr == as.numeric(substr(IL_filenames$data_chr, 4, 5))), ]$filename)
#     IL_id <- as.character(IL_filenames$IL[IL_filenames$IL == IL_list[i]][1])  ###### work on this part

    IL_df <- data.frame()
    for (j in 1:length(ILs_to_input)) {
        if (file.info(ILs_to_input[j])$size > 0) {
            IL_df_temp <- read.table(ILs_to_input[j])
            IL_df_temp <- cbind(IL_df_temp, unlist(strsplit(ILs_to_input[j], "\\.geno.*")))
            IL_df <- rbind(IL_df,IL_df_temp)
        }
    }
    


    if (nrow(IL_df) > 0) {
        colnames(IL_df) <- c("chr","pos","M", "P", "COV", "IL")

    IL_df$IL <- paste(sep="", "IL", gsub("\\.", "-", IL_df$IL))


        IL_df$COV.plot <- IL_df$COV #CHANGED
        IL_df$COV.plot[(IL_df$P-IL_df$M)/IL_df$COV < 0 & IL_df$COV > 0] <- -IL_df$COV.plot[(IL_df$P-IL_df$M)/IL_df$COV < 0 & IL_df$COV > 0] #CHANGED
        chr_formatted <- paste("chr", substr(IL_df$chr,9,10), sep = "")
        IL_df$chr <- chr_formatted
        max_COV <- max(IL_df$COV) #CHANGED        
        ylimits <- max_COV * offset
        geno <- ggplot(IL_df, aes(pos, COV.plot)) #CHANGED
        geno <- geno + geom_point(size = 1, aes(color = (P-M)/COV)) + scale_color_gradient2(low = "magenta", mid = "black", high = "green", name = expression((Pen - M82)/total)) + facet_grid(IL ~ .)
        # geno <- geno + opts(title = IL_df$chr[1])
        geno <- geno + guides(color = guide_legend(reverse=TRUE))
        geno <- geno + opts(title = paste( sep = " ", "Chromosome", as.numeric(substr(IL_df$chr[1],4,5)) ) )

        #added due to broken feature in ggplot 0.9.1:
        max_pos <- max(IL_df$pos)
        temp_max <- floor(max_pos/10000000)
        max_pos_fixed <- temp_max * 10000000
        label_max <- temp_max * 10
        geno <- geno + scale_x_continuous('Position on chromosome (Mb)', breaks = seq(0, max_pos_fixed, 10000000), labels = seq(0, label_max, 10)) + scale_y_continuous('Coverage', breaks = c(-100, 0, 100), labels = c("", "", ""), limits = c(-ylimits, ylimits)) #CHANGED

        geno <- geno + opts(panel.grid.minor=theme_blank(), panel.grid.major=theme_blank())

        geno <- geno + geom_segment(data = boundaries[boundaries$chr == chromosomes[i] & boundaries$IL %in% valid_ILs, ], facets = IL ~ ., color = "gray50", aes(x=start, xend=end, y = PEN * ylimits, yend = PEN * ylimits))


        # ggsave(paste("../plot_nr/split.ILs.", chr_formatted, ".pdf", sep = ""), plot = geno, width = 10.5, height = 8) #CHANGED
        ggsave(paste("RNASEQ_6x5.", chr_formatted, ".png", sep = ""), plot = geno, width = 6, height = 5) #CHANGED
        ggsave(paste("RNASEQ_6.8x4.8.", chr_formatted, ".png", sep = ""), plot = geno, width = 6.8, height = 4.8) #CHANGED
        # ggsave(paste("../plot/split.ILs.", chr_formatted, ".pdf", sep = ""), plot = geno, width = 10.5, height = 8) #CHANGED
        # ggsave(paste("../plot/split.ILs.", chr_formatted, ".png", sep = ""), plot = geno, width = 10.5, height = 8) #CHANGED
    }
}


#TODO
# customize IL order?
# sizes <- factor( c( "small", "large", "large", "small", "medium"),
#                  levels = c("small", "medium", "large"))





# > IL_filenames[(IL_filenames$IL_chr == 9) & (IL_filenames$data_chr == "chr12"),]
#                     filename    IL data_chr IL_chr
# 858 9.1.2.genotyped.chr12.nr 9.1.2    chr12      9
# 871 9.1.3.genotyped.chr12.nr 9.1.3    chr12      9
# 884   9.1.genotyped.chr12.nr   9.1    chr12      9
# 897 9.2.5.genotyped.chr12.nr 9.2.5    chr12      9
# 910 9.2.6.genotyped.chr12.nr 9.2.6    chr12      9
# 923   9.2.genotyped.chr12.nr   9.2    chr12      9
# 936 9.3.1.genotyped.chr12.nr 9.3.1    chr12      9
# 949 9.3.2.genotyped.chr12.nr 9.3.2    chr12      9
# 962   9.3.genotyped.chr12.nr   9.3    chr12      9

# > IL_filenames[936,]
#                     filename    IL data_chr IL_chr
# 936 9.3.1.genotyped.chr12.nr 9.3.1    chr12     9

# IL_filenames$IL_chr[936] = 12

# > IL_filenames[936,]
#                     filename    IL data_chr IL_chr
# 936 9.3.1.genotyped.chr12.nr 9.3.1    chr12     12


#then re-run w/ i = 4





