library("ggplot2")

offset <- 1.05

# added due to broken feature in ggplot 0.9.1:
# is it still broken/required?
max_pos <- max(geno_df$pos)
temp_max <- floor(max_pos / 10000000)
max_pos_fixed <- temp_max * 10000000
label_max <- temp_max * 10

max <-  max(geno_df$cov.plot)
min <- -min(geno_df$cov.plot)
if(max < 0) max <- 0
if(min < 0) min <- 0
max_lab <- (max %/% 10) * 10
min_lab <- (min %/% 10) * 10

geno.plot <- ggplot(geno_df, aes(pos, cov.plot)) +
  geom_point(size = 1, aes(color = (par2 - par1) / cov)) +
  facet_grid(chr ~ .) +
  scale_color_gradient2(
    low    = col_par1,
    mid    = col_het,
    high   = col_par2,
    limits = c(-1, 1),
    name   = substitute(frac(par2 - par1, total),
        list(par1 = par1, par2 = par2, total = "Total reads"))
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  ggtitle(id) +
  scale_x_continuous(
    'Position on chromosome (Mb)',
    breaks = seq(0, max_pos_fixed, 10000000),
    labels = seq(0, label_max,     10)
  ) +
  scale_y_continuous(
    'Coverage',
    breaks = c(-min_lab, 0, max_lab),
    labels = c( min_lab, 0, max_lab),
    limits = c(-offset*min, offset*max)
  )


