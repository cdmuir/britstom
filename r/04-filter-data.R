source("R/header.R")

stomata <- read_csv(str_c(path_proc_data, "/stomata.csv"))
phy <- read.nexus(file = str_c(path_proc_data, "/Lim_etal_2014_final.nex"))

##### Filter hydrophytes, helophytes, c4, cam, and nonangiosperm plants -----

nHydrophyte <- length(which(stomata$lifeform == "Hy"))
nC4 <- length(which(stomata$photo == "c4"))
nCAM <- length(which(stomata$photo == "cam"))

stomata %<>% filter(lifeform %in% c("Ph", "Ch", "hc", "Gn", "Th"),
                    photo == "c3") %>%
  select(-photo)

angio_phy <- extract.clade(phy, node = length(phy$tip.label) + 
                             which(phy$node.label == "angiosperms"))

stomata %<>% filter(species %in% angio_phy$tip.label)

# Export new dataset and tree
write_csv(stomata, str_c(path_proc_data, "/stomata_filtered.csv"))
write.nexus(angio_phy, file = str_c(path_proc_data, "/angio_phy.nex"))
