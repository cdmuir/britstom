source("R/header.R")

phy <- read.nexus(file = str_c(path_proc_data, "/angio_phy_modified.nex"))
stomata <- read_csv(str_c(path_proc_data, "/stomata_filtered.csv"))

# Arrange stomata in same order as phylogeny for plotting
stomata %<>% extract(match(phy$tip.label, .$species), 1:ncol(.))

lf <- c("chamaephyte", "geophyte", "hemicryptophyte", "hydrophyte", "phanerophyte", "therophyte")
names(lf) <- c("Ch", "Gn", "hc", "Hy", "Ph", "Th")
lf <- lf[names(sort(tapply(stomata$sr_propAd, stomata$lifeform, mean)))]

# Factor lifeform for plotting
stomata$lifeform %<>% factor(levels = names(lf))

# Rescale phylogeny for plotting
phy$edge.length %<>% divide_by(max(nodeHeights(phy)[, 2]))

# Angle for calculating position of data in fan phylogeny
theta_mid <- seq(0, 2 * pi, length.out = Ntip(phy) + 1)
theta_mid %<>% extract(1:Ntip(phy))
theta_left <- theta_mid - pi / Ntip(phy)
theta_right <- theta_mid + pi / Ntip(phy)

# Main plot
pdf(str_c(path_figures, "/figure_phylo.pdf"), w = 6.5, h = 6.5)#, 
   # useDingbats = FALSE)
par(mai = c(1, 0, 1, 2))
gp <- plot(phy, type = "fan", show.tip.label = FALSE,
           x.lim = c(-1.25, 1.25), y.lim = c(-1.25, 1.25))

# Growth form

  ## Palette
  palette(brewer.pal(5, "Greens"))

  ## Plot wedges
  message("double check that wedges are drawn in same order as data")
  for (i in 1:Ntip(phy)) {
    draw_wedge(theta_left[i], theta_right[i], r1 = 1.025, r2 = 1.075, 
               col = stomata$lifeform[i])
  }  

  ## Legend
  clip(grconvertX(0, "ndc", "user"), grconvertX(1, "ndc", "user"),
       grconvertY(0, "ndc", "user"), grconvertY(1, "ndc", "user"))
  text(grconvertX(5.5, "in", "user"), grconvertY(0.85, "ndc", "user"),
       labels = "Growth form", cex = 1.5)
  points(rep(grconvertX(5, "in", "user"), 5),
         seq(0.1, 0.9, 0.2) * 1.25,
         pch = 21, col = "black", bg = 1:5, cex = 3)
  text(rep(grconvertX(5, "in", "user"), 5),
       seq(0.1, 0.9, 0.2) * 1.25,
       labels = lf, pos = 4, offset = 1)

# Stomatal ratio
  
  ## Histofan
  draw_histofan(stomata$sr_even, 1.1, 1.2)
  
  ## Legend
  text(grconvertX(5.5, "in", "user"), grconvertY(0.45, "ndc", "user"),
       labels = "Stomatal ratio", cex = 1.5)
  
  theta <- seq(-pi / 12, pi / 12, length.out = 1e2)
  r <- grconvertX(6, "in", "user")
  x <- cos(theta) * r
  y <- sin(theta) * r + grconvertY(0.25, "ndc", "user")
  points(x, y, type = "l", col = "grey")
  text(x[1], y[1], labels = "1", pos = 1, 
       srt = 360 * theta[1] / (2 * pi))
  
  r <- grconvertX(5.5, "in", "user")
  x <- cos(theta) * r
  y <- sin(theta) * r + grconvertY(0.25, "ndc", "user")
  points(x, y, type = "l", col = "grey")
  text(x[1], y[1], labels = "0.5", pos = 1, 
       srt = 360 * theta[1] / (2 * pi))
  text(x[1], y[1], labels = expression(SR[even]), pos = 1, offset = 2, 
       cex = 1.25, srt = 360 * theta[1] / (2 * pi))
  
  r <- grconvertX(5, "in", "user")
  x <- cos(theta) * r
  y <- sin(theta) * r + grconvertY(0.25, "ndc", "user")
  points(x, y, type = "l", col = "grey")
  text(x[1], y[1], labels = "0", pos = 1, 
       srt = 360 * theta[1] / (2 * pi))
  
  set.seed(883090777)
  d <- runif(8, grconvertX(5, "in", "user"), grconvertX(6, "in", "user"))
  theta <- seq(from = min(theta), to = max(theta), length.out = length(d) + 1)
  X <- Y <- numeric()
  for (i in 1:length(d)) {
    X %<>% c(cos(seq(theta[i], theta[i + 1], length.out = 1e2)) * d[i])
    Y %<>% c(sin(seq(theta[i], theta[i + 1], length.out = 1e2)) * d[i])
  }
  
  Y %<>% add(grconvertY(0.25, "ndc", "user"))
  points(X, Y, type = "l", lwd = 3)
  
  dev.off()
  
