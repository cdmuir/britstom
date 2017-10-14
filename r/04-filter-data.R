source("r/header.R")

stomata <- read_csv(str_c(path_proc_data, "/stomata.csv"), 
                    col_types = cols(
                      species = col_character(),
                      ab_density = col_double(),
                      ad_density = col_double(),
                      acceptedname = col_character(),
                      sr_propAd = col_double(),
                      sr_even = col_double(),
                      photo = col_character(),
                      lifeform = col_character(),
                      ellenberg_light = col_integer(),
                      habit = col_character()
                    ))
phy <- read.nexus(file = str_c(path_proc_data, "/Lim_etal_2014_final.nex"))

##### Filter hydrophytes, helophytes, c4, cam, and nonangiosperm plants -----

nHydrophyte <- length(which(stomata$lifeform == "Hy"))
nC4 <- length(which(stomata$photo == "c4"))
nCAM <- length(which(stomata$photo == "cam"))

stomata %<>% filter(lifeform %in% c("Ph", "Ch", "hc", "Gn", "Th"),
                    photo == "c3") %>%
  dplyr::select(-photo)

angio_phy <- extract.clade(phy, node = length(phy$tip.label) + 
                             which(phy$node.label == "angiosperms"))

# Change species to acceptedname is latter is in phylogeny
stomata$species[which(!stomata$species %in% phy$tip.label)] <-
  stomata$acceptedname[which(!stomata$species %in% phy$tip.label)]
stomata %<>% dplyr::select(-acceptedname)

n <- nrow(stomata)

stomata %<>% filter(species %in% angio_phy$tip.label)

nNonangio <- n - nrow(stomata)
rm(n)

# Export new dataset and tree
write_csv(stomata, str_c(path_proc_data, "/stomata_filtered.csv"))
write.nexus(angio_phy, file = str_c(path_proc_data, "/angio_phy.nex"))

# Export objects to ms
export2ms(c("nC4", "nCAM", "nHydrophyte", "nNonangio", "stomata"))
