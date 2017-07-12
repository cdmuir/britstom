source("R/header.R")

phy <- read.nexus(file = str_c(path_proc_data, "/angio_phy_modified.nex"))
stomata <- read_csv(str_c(path_proc_data, "/stomata_filtered.csv"))

##### Phylogenetic structural equation model ---
  
# small data.frame for Rphylopars
td <- stomata %>% select(species, ab_density, ad_density, ellenberg_light) %>%
  mutate(logSDab = log(ab_density + 1), logSDad = log(ad_density + 1)) %>%
  select(-ab_density, -ad_density) %>%
  as.data.frame() %>% 
  set_rownames(.$species) # need to change row names for phylolm

# ppOU <- phylopars(td, phy, model = "mvOU")
# write_rds(ppOU, str_c(path_objects, "/ppOU.rds"))
ppOU <- read_rds(str_c(path_objects, "/ppOU.rds"))
phy_mat <- ppOU$pars$phylocov
  
# define (co)variances for use in text
varSDab <- phy_mat["logSDab", "logSDab"]
varSDad <- phy_mat["logSDad", "logSDad"]
covSD <- phy_mat["logSDab", "logSDad"]
varSR <- varSDab + varSDad - 2 * covSD
  
ppBM <- phylopars(td, phy, model = "BM")
  
# AIC comparison of OU and BM
AIC(ppOU, ppBM)
  
# Strucutural equation model using lavaan
model <- 'logSDab ~ ellenberg_light
          logSDad ~ ellenberg_light'
  
fitSEM <- sem(model, sample.cov = ppOU$pars$phylocov, sample.nobs = nrow(td))
  
# Calculate p-values based on z-scores
pvalSEM <- pnorm(fitSEM@ParTable$est / fitSEM@ParTable$se, lower.tail = FALSE, log = TRUE)
names(pvalSEM) <- paste(fitSEM@ParTable$lhs, fitSEM@ParTable$op, fitSEM@ParTable$rhs)
