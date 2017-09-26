source("R/header.R")

phy <- read.nexus(file = str_c(path_proc_data, "/angio_phy_modified.nex"))
stomata <- read_csv(str_c(path_proc_data, "/stomata_filtered.csv"))

# Arrange stomata in same order as phylogeny for plotting
stomata %<>% extract(match(phy$tip.label, .$species), 1:ncol(.))

# Add blank rows for internal node data
blank_rows <- matrix(NA, nrow = nrow(stomata) - 1, ncol = ncol(stomata)) %>%
  as.data.frame() %>%
  set_colnames(colnames(stomata))

stomata %<>% bind_rows(blank_rows)
rm(blank_rows)

##### Raunikaer life form and Ellenberg light indicator values versus sr_even -----

lf <- c("chamaephyte", "geophyte", "hemicryptophyte", "hydrophyte", "phanerophyte", "therophyte")
names(lf) <- c("Ch", "Gn", "hc", "Hy", "Ph", "Th")
lf <- lf[names(sort(tapply(stomata$sr_propAd, stomata$lifeform, mean)))]

plot(phy, type = "fan", show.tip.label = FALSE)


##### Scratch -----

stomata$lifeform %<>% 
  factor(levels = names(lf))
gp <- ggtree(phy, layout = "circular") +
  geom_tippoint(aes(color = stomata$lifeform), shape = 18, size = 3) +
  scale_colour_brewer("Growth form", type= "seq", palette = "Greens") +
  theme(legend.position = "right")
print(gp)
