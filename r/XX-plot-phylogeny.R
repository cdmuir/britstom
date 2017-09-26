source("R/header.R")

phy <- read.nexus(file = str_c(path_proc_data, "/angio_phy_modified.nex"))
stomata <- read_csv(str_c(path_proc_data, "/stomata_filtered.csv"))

# Arrange stomata in same order as phylogeny for plotting
stomata %<>% extract(match(phy$tip.label, .$species), 1:ncol(.))

# Add blank rows for internal node data (this is for ggtree - remove if unused)
# blank_rows <- matrix(NA, nrow = nrow(stomata) - 1, ncol = ncol(stomata)) %>%
#   as.data.frame() %>%
#   set_colnames(colnames(stomata))
# 
# stomata %<>% bind_rows(blank_rows)
# rm(blank_rows)

##### Raunikaer life form and Ellenberg light indicator values versus sr_even -----

lf <- c("chamaephyte", "geophyte", "hemicryptophyte", "hydrophyte", "phanerophyte", "therophyte")
names(lf) <- c("Ch", "Gn", "hc", "Hy", "Ph", "Th")
lf <- lf[names(sort(tapply(stomata$sr_propAd, stomata$lifeform, mean)))]

# Factor lifeform for plotting
stomata$lifeform %<>% factor(levels = names(lf))

##### Scratch -----

# gp <- ggtree(phy, layout = "rectangular") +
#   geom_tippoint(aes(color = stomata$lifeform), shape = 19, size = 3) +
#   scale_colour_brewer("Growth form", type= "seq", palette = "Greens") +
#   theme(legend.position = "right")
# print(gp)


#### I'm liking this appraoch. ggtree not flexible enough ####

# Rescale phylogeny for plotting
phy$edge.length %<>% divide_by(max(nodeHeights(phy)[, 2]))

# Angle for calculating position of data in fan phylogeny
theta_mid <- seq(0, 2 * pi, length.out = Ntip(phy))
theta_left <- theta_mid - pi / Ntip(phy)
theta_right <- theta_mid + pi / Ntip(phy)

# Main plot
pdf(str_c(path_figures, "/figure_phylo.pdf"), w = 6.5, h = 6.5, 
    useDingbats = FALSE)
par(mai = c(1, 0, 1, 2))
gp <- plot(phy, type = "fan", show.tip.label = FALSE,
           x.lim = c(-1.25, 1.25), y.lim = c(-1.25, 1.25))

# Growth form

  ## Palette
  palette(brewer.pal(5, "Greens"))

  ## Plot wedges
  for (i in 1:Ntip(phy)) {
    draw_wedge(theta_left[i], theta_right[i], r1 = 1.025, r2 = 1.075, 
               col = stomata$lifeform[i])
  }  

  ## Legend? ...
  dev.off()
  

  x_gf <- cos(theta) * r_gf
  y_gf <- sin(theta) * r_gf
  
  ## Plot points
  points(x_gf, y_gf, col = stomata$lifeform)

  y <- sin(theta) * 1.02
x <- cos(theta) * 1.02
points(x, y, pch = ".")

# legend - get glyphs for each growth form