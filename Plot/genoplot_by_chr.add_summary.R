boundaries <- read.table("IL_boundaries", head = T)

# > head (boundaries)
#   chr      IL    start      end par2
# 1   1   IL1-1        1   303076  -1
# 2   1   IL1-1   309021 77720421   1
# 3   1   IL1-1 77721171 90304244  -1
# 4   1 IL1-1-2        1   303076  -1
# 5   1 IL1-1-2   309061  1137671   1
# 6   1 IL1-1-2  1151331 90304244  -1

geno.plot <- geno.plot + geom_segment(
    data   = boundaries[boundaries$chr == chromosomes[i] & boundaries$IL %in% valid_ILs, ],
    facets = IL ~ .,
    color  = "gray50",
    aes(
      x    = start,
      xend = end,
      y    = par2 * ylimits,
      yend = par2 * ylimits)
  )