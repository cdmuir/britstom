source("r/header.R")

stomata <- read_csv(str_c(path_proc_data, "/stomata_filtered.csv"), 
                    col_types = cols(
                      species = col_character(),
                      ab_density = col_double(),
                      ad_density = col_double(),
                      sr_propAd = col_double(),
                      sr_even = col_double(),
                      lifeform = col_character(),
                      ellenberg_light = col_integer(),
                      growthform = col_character()
                    ))

gf <- c("tree", "shrub", "perennial", "biennial", "annual")

lf <- c("chamaephyte", "geophyte", "hemicryptophyte", "hydrophyte", "phanerophyte", "therophyte")
names(lf) <- c("Ch", "Gn", "hc", "Hy", "Ph", "Th")
lf <- lf[names(sort(tapply(stomata$sr_propAd, stomata$lifeform, mean)))]

palette(brewer.pal(5, "Greens"))

palette(brewer.pal(5, "Greens"))

pdf(str_c(path_figures, "/figureS_compare-growth-forms.pdf"), 
    width = 4, height = 6, useDingbats = FALSE)

  par(mar = c(5, 5, 1, 1))
  tab <- prop_table(stomata$lifeform, stomata$growthform)
  tab <- tab[, names(lf)]
  tab %<>% set_colnames(lf)
  dotchart(tab, bg = 1:5, pt.cex = 1.5, xlab = "Proportion", cex.lab = 1.5)

dev.off()
