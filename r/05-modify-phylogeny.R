source("R/header.R")

stomata <- read_csv(str_c(path_proc_data, "/stomata_filtered.csv"))
angio_phy <- read.nexus(file = str_c(path_proc_data, "/angio_phy.nex"))

##### Make dichotomous and slightly adjust edge lengths to make exactly ultrametric -----

angio_phy %<>% multi2di()
angio_phy$edge.length[angio_phy$edge.length == 0] <- 0.02
x <- diag(vcv.phylo(angio_phy))
is.ultrametric(angio_phy)
n <- length(angio_phy$tip.label)
te <- sapply(1:n, function(x, y) which(y == x), y = angio_phy$edge[, 2]) # terminal edges
angio_phy$edge.length[te] <- angio_phy$edge.length[te] + (max(x) - x)
is.ultrametric(angio_phy)
angio_phy %<>% drop.tip(.$tip.label[!.$tip.label %in% stomata$species])

write.nexus(angio_phy, file = str_c(path_proc_data, "/angio_phy_modified.nex"))
