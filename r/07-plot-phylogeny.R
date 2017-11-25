source("R/header.R")

phy <- read.nexus(file = str_c(path_proc_data, "/angio_phy_modified.nex"))
stomata <- read_csv(str_c(path_proc_data, "/stomata_filtered.csv"), 
                    col_types = cols(
                      species = col_character(),
                      ab_density = col_double(),
                      ad_density = col_double(),
                      sr_propAd = col_double(),
                      sr_even = col_double(),
                      lifeform = col_character(),
                      ellenberg_light = col_integer(),
                      habit = col_character()
                    ))

# Arrange stomata in same order as phylogeny for plotting
stomata %<>% magrittr::extract(match(phy$tip.label, .$species), 1:ncol(.))

hf <- c("tree", "shrub", "perennial", "biennial", "annual")

lf <- c("chamaephyte", "geophyte", "hemicryptophyte", "phanerophyte", "therophyte")
names(lf) <- c("Ch", "Gn", "hc", "Ph", "Th")
lf <- lf[names(sort(tapply(stomata$sr_propAd, stomata$lifeform, mean)))]

# Factor habit and lifeform for plotting
stomata$habit %<>% factor(levels = hf)
stomata$lifeform %<>% factor(levels = names(lf))

# Rescale phylogeny for plotting
phy$edge.length %<>% divide_by(max(nodeHeights(phy)[, 2]))

# Angle for calculating position of data in fan phylogeny
theta_mid <- seq(0, 2 * pi, length.out = Ntip(phy) + 1)
theta_mid %<>% magrittr::extract(1:Ntip(phy))
theta_left <- theta_mid - pi / Ntip(phy)
theta_right <- theta_mid + pi / Ntip(phy)

##### Using habit -----

# Main plot
pdf(str_c(path_figures, "/figure_phylo-habit.pdf"), w = 6.5, h = 6.5, 
    useDingbats = FALSE)
par(mai = c(1, 0, 1, 2))
gp <- plot(phy, type = "fan", show.tip.label = FALSE,
           x.lim = c(-1.25, 1.25), y.lim = c(-1.25, 1.25))

# Growth form

  ## Palette for growth form
  palette(brewer.pal(5, "Greens"))
  palette(rev(palette()))
  
  ## Plot wedges
  for (i in 1:Ntip(phy)) {
    draw_wedge(theta_left[i], theta_right[i], r1 = 1.025, r2 = 1.075, 
               col = stomata$habit[i])
  }  

  ## Legend
  clip(grconvertX(0, "ndc", "user"), grconvertX(1, "ndc", "user"),
       grconvertY(0, "ndc", "user"), grconvertY(1, "ndc", "user"))
  text(grconvertX(5.5, "in", "user"), grconvertY(0.9, "ndc", "user"),
       labels = "Growth form", cex = 1.5)
  points(rep(grconvertX(5, "in", "user"), 5),
         seq(0.15, 0.95, 0.2) * grconvertY(0.85, "ndc", "user"),
         pch = 21, col = "black", bg = 1:5, cex = 3)
  text(rep(grconvertX(5, "in", "user"), 5),
       seq(0.15, 0.95, 0.2) * grconvertY(0.85, "ndc", "user"),
       labels = hf, pos = 4, offset = 1)
  
# Light
  
  ## Palette for light
  palette(brewer.pal(length(unique(stomata$ellenberg_light)), "Oranges"))
  palette(rev(palette()))
  
  ## Plot wedges
  for (i in 1:Ntip(phy)) {
    draw_wedge(theta_left[i], theta_right[i], r1 = 1.075, r2 = 1.125, 
               col = stomata$ellenberg_light[i] - min(stomata$ellenberg_light) + 1)
  }  
  
  ## Legend
  clip(grconvertX(0, "ndc", "user"), grconvertX(1, "ndc", "user"),
       grconvertY(0, "ndc", "user"), grconvertY(1, "ndc", "user"))
  text(grconvertX(5.5, "in", "user"), grconvertY(0.45, "ndc", "user"),
       labels = "Ellenberg light\nindicator value", cex = 1.5)
  points(rep(grconvertX(5, "in", "user"), 7),
         seq(-0.4, -1.2, length.out = 7) * 1.25,
         pch = 21, col = "black", bg = 1:7, cex = 3)
  text(rep(grconvertX(5, "in", "user"), 5),
       seq(-0.4, -1.2, length.out = 7) * 1.25,
       labels = str_c(3:9, c(" shade", rep("", 5), " sun")), pos = 4, offset = 1)
  
# Stomatal ratio

  ## Histofan
  draw_histofan(stomata$sr_even, 1.15, 1.3)

  ## Legend
  text(grconvertX(0.1, "ndc", "user"), grconvertY(0.9, "ndc", "user"),
       labels = "Stomatal\nratio", cex = 1.5)

  theta <- seq(1 / 2 * pi - pi / 12, 1 / 2 * pi + pi / 12, length.out = 1e2)
  theta_axis <- theta[1] - 24 / 360
  theta_lab <- theta[1] - 72 / 360
  
  r <- grconvertX(5.25, "in", "user")
  x <- cos(theta) * r
  y <- sin(theta) * r
  points(x, y, type = "l", col = "grey")
  text(cos(theta_axis) * r, sin(theta_axis) * r, 
       labels = "1", srt = 360 * (pi - theta_axis / (2 * pi)))
  
  r <- grconvertX(5, "in", "user")
  x <- cos(theta) * r
  y <- sin(theta) * r
  points(x, y, type = "l", col = "grey")
  text(cos(theta_axis) * r, sin(theta_axis) * r, 
       labels = "0.5", srt = 360 * (pi - theta_axis / (2 * pi)))
  text(cos(theta_lab) * r, sin(theta_lab) * r, 
       labels = expression(SR[even]), cex = 1.25,
       srt = 360 * (theta_lab / (2 * pi)))

  r <- grconvertX(4.75, "in", "user")
  x <- cos(theta) * r
  y <- sin(theta) * r
  points(x, y, type = "l", col = "grey")
  text(cos(theta_axis) * r, sin(theta_axis) * r,
       labels = "0", srt = 360 * (pi - theta_lab / (2 * pi)))

  set.seed(883090777)
  d <- runif(8, grconvertX(4.75, "in", "user"), grconvertX(5.25, "in", "user"))
  theta <- seq(from = min(theta), to = max(theta), length.out = length(d) + 1)
  X <- Y <- numeric()
  for (i in 1:length(d)) {
    X %<>% c(cos(seq(theta[i], theta[i + 1], length.out = 1e2)) * d[i])
    Y %<>% c(sin(seq(theta[i], theta[i + 1], length.out = 1e2)) * d[i])
  }

  points(X, Y, type = "l", lwd = 3)

dev.off()

##### Using Raunkiaer lifeform -----

# Main plot
pdf(str_c(path_figures, "/figure_phylo-lifeform.pdf"), w = 6.5, h = 6.5, 
    useDingbats = FALSE)
par(mai = c(1, 0, 1, 2))
gp <- plot(phy, type = "fan", show.tip.label = FALSE,
           x.lim = c(-1.25, 1.25), y.lim = c(-1.25, 1.25))

# Growth form

  ## Palette
  palette(brewer.pal(5, "Greens"))
  palette(rev(palette()))

  ## Plot wedges
  for (i in 1:Ntip(phy)) {
    draw_wedge(theta_left[i], theta_right[i], r1 = 1.025, r2 = 1.075, 
               col = stomata$lifeform[i])
  }  

  ## Legend
  clip(grconvertX(0, "ndc", "user"), grconvertX(1, "ndc", "user"),
       grconvertY(0, "ndc", "user"), grconvertY(1, "ndc", "user"))
  text(grconvertX(5.5, "in", "user"), grconvertY(0.9, "ndc", "user"),
       labels = "Growth form", cex = 1.5)
  points(rep(grconvertX(5, "in", "user"), 5),
         seq(0.15, 0.95, 0.2) * grconvertY(0.85, "ndc", "user"),
         pch = 21, col = "black", bg = 1:5, cex = 3)
  text(rep(grconvertX(5, "in", "user"), 5),
       seq(0.15, 0.95, 0.2) * grconvertY(0.85, "ndc", "user"),
       labels = lf, pos = 4, offset = 1)

# Light
  
  ## Palette for light
  palette(brewer.pal(length(unique(stomata$ellenberg_light)), "Oranges"))
  palette(rev(palette()))
  
  ## Plot wedges
  for (i in 1:Ntip(phy)) {
    draw_wedge(theta_left[i], theta_right[i], r1 = 1.075, r2 = 1.125, 
               col = stomata$ellenberg_light[i] - min(stomata$ellenberg_light) + 1)
  }  
  
  ## Legend
  clip(grconvertX(0, "ndc", "user"), grconvertX(1, "ndc", "user"),
       grconvertY(0, "ndc", "user"), grconvertY(1, "ndc", "user"))
  text(grconvertX(5.5, "in", "user"), grconvertY(0.45, "ndc", "user"),
       labels = "Ellenberg light\nindicator value", cex = 1.5)
  points(rep(grconvertX(5, "in", "user"), 7),
         seq(-0.4, -1.2, length.out = 7) * 1.25,
         pch = 21, col = "black", bg = 1:7, cex = 3)
  text(rep(grconvertX(5, "in", "user"), 5),
       seq(-0.4, -1.2, length.out = 7) * 1.25,
       labels = str_c(3:9, c(" shade", rep("", 5), " sun")), pos = 4, offset = 1)
  
# Stomatal ratio
  
  ## Histofan
  draw_histofan(stomata$sr_even, 1.15, 1.3)
  
  ## Legend
  text(grconvertX(0.1, "ndc", "user"), grconvertY(0.9, "ndc", "user"),
       labels = "Stomatal\nratio", cex = 1.5)
  
  theta <- seq(1 / 2 * pi - pi / 12, 1 / 2 * pi + pi / 12, length.out = 1e2)
  theta_axis <- theta[1] - 24 / 360
  theta_lab <- theta[1] - 72 / 360
  
  r <- grconvertX(5.25, "in", "user")
  x <- cos(theta) * r
  y <- sin(theta) * r #+ grconvertY(0.25, "ndc", "user")
  points(x, y, type = "l", col = "grey")
  text(cos(theta_axis) * r, sin(theta_axis) * r, 
       labels = "1", srt = 360 * (pi - theta_axis / (2 * pi)))
  
  r <- grconvertX(5, "in", "user")
  x <- cos(theta) * r
  y <- sin(theta) * r #+ grconvertY(0.25, "ndc", "user")
  points(x, y, type = "l", col = "grey")
  text(cos(theta_axis) * r, sin(theta_axis) * r, 
       labels = "0.5", srt = 360 * (pi - theta_axis / (2 * pi)))
  text(cos(theta_lab) * r, sin(theta_lab) * r, 
       labels = expression(SR[even]), cex = 1.25,
       srt = 360 * (theta_lab / (2 * pi)))
  
  r <- grconvertX(4.75, "in", "user")
  x <- cos(theta) * r
  y <- sin(theta) * r #+ grconvertY(0.25, "ndc", "user")
  points(x, y, type = "l", col = "grey")
  text(cos(theta_axis) * r, sin(theta_axis) * r,
       labels = "0", srt = 360 * (pi - theta_lab / (2 * pi)))
  
  set.seed(883090777)
  d <- runif(8, grconvertX(4.75, "in", "user"), grconvertX(5.25, "in", "user"))
  theta <- seq(from = min(theta), to = max(theta), length.out = length(d) + 1)
  X <- Y <- numeric()
  for (i in 1:length(d)) {
    X %<>% c(cos(seq(theta[i], theta[i + 1], length.out = 1e2)) * d[i])
    Y %<>% c(sin(seq(theta[i], theta[i + 1], length.out = 1e2)) * d[i])
  }
  
  points(X, Y, type = "l", lwd = 3)
  
dev.off()
  
