source("R/header.R")

phy <- read.nexus(file = str_c(path_proc_data, "/angio_phy_modified.nex"))
stomata <- read_csv(str_c(path_proc_data, "/stomata_filtered.csv"))

##### Phylogenetic structural equation model ---
  
  # small data.frame for Rphylopars
td <- stomata

  td <- data.frame(species = stomata$species)
  td$logSDab <- log(stomata$ab_density + 1)
  td$logSDad <- log(stomata$ad_density + 1)
  td$ellenberg_light <- stomata$ellenberg_light
  
  # ppOU <- phylopars(td, modAngioPhy, model = "mvOU")
  # saveRDS(ppOU, paste0(pathR, "/ppOU.rds"))
  ppOU <- readRDS(paste0(pathR, "/ppOU.rds"))
  phyMat <- ppOU$pars$phylocov
  
  # define (co)variances for use in text
  varSDab <- phyMat["logSDab", "logSDab"]
  varSDad <- phyMat["logSDad", "logSDad"]
  covSD <- phyMat["logSDab", "logSDad"]
  varSR <- varSDab + varSDad - 2 * covSD
  
  ppBM <- phylopars(td, modAngioPhy, model = "BM")
  
  # AIC comparison of OU and BM
  AIC(ppOU, ppBM)
  
  # Strucutural equation model using lavaan
  model <- 'logSDab ~ ellenberg_light
  logSDad ~ ellenberg_light'
  
  fitSEM <- sem(model, sample.cov = ppOU$pars$phylocov, sample.nobs = nrow(td))
  
  # Calculate p-values based on z-scores
  pvalSEM <- pnorm(fitSEM@ParTable$est / fitSEM@ParTable$se, lower.tail = FALSE, log = TRUE)
  names(pvalSEM) <- paste(fitSEM@ParTable$lhs, fitSEM@ParTable$op, fitSEM@ParTable$rhs)
  
  rm("td")
  
  @