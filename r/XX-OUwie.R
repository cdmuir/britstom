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

library(OUwie)

# Import fitted model
fitSR_lf_pl <- read_rds(str_c(path_objects, "/fitSR_lf_pl.rds"))

# Prepare phylogeny and data for OUwie
phy$node.label <- rep("1", Nnode(phy))
OUwie_dat <- data.frame(species = phy$tip.label, regime = "1", X = 0)

#OUwie.sim(phy, alpha = fitSR_lf_a$optpar, sigma.sq = fitSR_lf_a$sigma2,
#          theta = c(0, 0), theta0 = 0)

# Sample root state randomly from stationary distribution
nsim <- 1e1
rnorm(nsim, 0, sqrt(fitSR_lf_pl$sigma2 / (2 * fitSR_lf_pl$optpar)))
sim_data <- OUwie.sim(phy, data = OUwie_dat, alpha = rep(fitSR_lf_pl$optpar, 2), 
                      sigma.sq = rep(fitSR_lf_pl$sigma2, 2),
                      theta = c(0, 0), theta0 = 0)
sd(sim_data$X)
sd(resid(fitSR_lf_pl))
