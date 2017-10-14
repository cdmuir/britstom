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

# Import fitted model
fitSR_lf_pl <- read_rds(str_c(path_objects, "/fitSR_lf_pl.rds"))

# Prepare phylogeny and data for OUwie
phy$node.label <- rep("1", Nnode(phy))
OUwie_dat <- data.frame(species = phy$tip.label, regime = "1", x = 0)

# Parameteric bootstrap parameters
par <- fitSR_lf_pl$bootstrap[, c("sigma2", "alpha")]

# Simulate data (takes several minutes)
# set.seed(493151857)
# actual_tree <- apply(par, 1, function(x) {
#   tmp <- OUwie.sim(phy, data = OUwie_dat, alpha = rep(x["alpha"], 2), 
#                    sigma.sq = rep(x["sigma2"], 2), theta = c(0, 0), theta0 = 0)
#   var(tmp$X)
# })
# 
# stationary_tree <- apply(par, 1, function(x) {
#   var(rnorm(Ntip(phy), 0, sqrt(x["sigma2"] / (2 * x["alpha"]))))
# })
# 
# var_comp <- data.frame(actual_tree, stationary_tree)
# write_rds(var_comp, str_c(path_objects, "/var_comp.rds"))

var_comp <- read_rds(str_c(path_objects, "/var_comp.rds"))

ttest <- t.test(var_comp$actual_tree, var_comp$stationary_tree, paired = TRUE)

# Export objects to ms
export2ms(c("ttest", "var_comp"))
