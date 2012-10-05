geno.plot <- geno.plot + geom_segment(
    data   = boundaries[boundaries$chr == chromosomes[i] & boundaries$IL %in% valid_ILs, ],
    facets = IL ~ .,
    color  = "gray50",
    aes(
      x    = start,
      xend = end,
      y    = PEN * ylimits,
      yend = PEN * ylimits)
  )