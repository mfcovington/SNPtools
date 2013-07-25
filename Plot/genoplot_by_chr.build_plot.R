# TODO:
# - add color customization
# - add plot size customization
# - add plot format customization

library("ggplot2")

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
  )