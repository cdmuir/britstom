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

##### Import phylogeny of British Flora -----
# from Lim et al 2014
tmp <- read_lines(str_c(path_raw_data, "/S15105.nex"))
l1 <- str_detect(tmp, "BEGIN CHARACTERS") %>% which() - 1
l2 <- str_detect(tmp, "BEGIN TREES") %>% which() - 1

# Lithospermum purpurocaeruleum is mispelled
tmp %<>% str_replace_all("Lithospermum_purpureocaeruleum", "Lithospermum_purpurocaeruleum")

# Changing Populus_nigra_sens.lat. to Populus nigra officinale 
tmp %<>% str_replace_all("Populus_nigra_sens.lat.", "Populus_nigra")

# Changing Rosa_canina_agg. to Rosa_canina officinale 
tmp %<>% str_replace_all("Rosa_canina_agg.", "Rosa_canina")

write_lines(tmp[c(1:l1, l2:length(tmp))], 
            path = str_c(path_raw_data, "/Lim_etal_2014.nex"))
phy <- read_nexus_phylo(str_c(path_raw_data, "/Lim_etal_2014.nex"))

##### Check for species missing from phylogeny -----
notInPhy <- which(!(stomata$species %in% phy$tip.label |
                      stomata$acceptedname %in% phy$tip.label))
stomata[notInPhy, ]
nNotInPhy1 <- length(notInPhy)

##### Modify tip labels -----

# Using Cerastium brachypetalum as proxy for C. cerastoides based on similar divergence in Scheen et al. 2004 (Am J Bot, 91(6): 943-952)
tmp %<>% str_replace_all("Cerastium_brachypetalum", "Cerastium_cerastoides")

# Using Coincya monensis as proxy for C. wrightii 
tmp %<>% str_replace_all("Coincya_monensis", "Coincya_wrightii")

# Using Gnaphalium luteoalbum as proxy for G. norvegicum 
tmp %<>% str_replace_all("Gnaphalium_luteoalbum", "Gnaphalium_norvegicum")

# Using Lythrum portula as proxy for L. hyssopifolia 
tmp %<>% str_replace_all("Lythrum_portula", "Lythrum_hyssopifolia")

# Using Peucedanum ostruthium as proxy for P. officinale 
tmp %<>% str_replace_all("Peucedanum_ostruthium", "Peucedanum_officinale")

# Using Stachys alpina as proxy for S. germanica based on similar divergence in Salmaki et al. 2013 (Mol Phylo Evol, 69(3): 535-551)
tmp %<>% str_replace_all("Stachys_alpina", "Stachys_germanica")

# Using Tripleurospermum_maritimum as proxy for T. inodorum 
tmp %<>% str_replace("Tripleurospermum_maritimum", "Tripleurospermum_inodorum")

# Using Viola_arvensis as proxy for V. kitaibeliana 
tmp %<>% str_replace_all("Viola_arvensis", "Viola_kitaibeliana")

# Write first modified phylogeny
write_lines(tmp[c(1:l1, l2:length(tmp))], 
            path = str_c(path_proc_data, "/Lim_etal_2014_mod1.nex"))
phy <- read_nexus_phylo(paste0(path_proc_data, "/Lim_etal_2014_mod1.nex"))

##### Bind additional tips -----
tips <- phy$tip.label
nodes <- sapply(tips, function(x,y) which(y == x), y = phy$tip.label)
edge.lengths <- setNames(phy$edge.length[sapply(nodes, function(x, y) {
  which(y == x) }, y = phy$edge[, 2])], names(nodes))

##### Check for species missing from phylogeny -----
notInPhy <- which(!(stomata$species %in% phy$tip.label |
                      stomata$acceptedname %in% phy$tip.label))
stomata[notInPhy, ]
nNotInPhy2 <- length(notInPhy)

# Make Agrostemma githago sister to all (Silene sp., Lychnis sp.) [based on Fior et al. 2006, Am J Bot 104(2):399-411]
n <- getMRCA(phy, phy$tip.label[c(grep("Silene", phy$tip.label), grep("Lychnis", phy$tip.label))])

# Branch length is average of distance from (Silene sp., Lychnis sp.) MRCA to (Silene sp., Lychnis sp.) sp.
el1 <- sapply(nodes[c(grep("Silene_", phy$tip.label), grep("Lychnis_", phy$tip.label))], 
              nodeheight, tree = phy)
el2 <- nodeheight(tree = phy, node = n)

phy <- bind.tip(phy, "Agrostemma_githago", edge.length = mean(el1 - el2), where = n)

# Make Heracleum mantegazzianum as sister to H. sphondylium
phy <- bind.tip(phy, "Heracleum_mantegazzianum", 
                edge.length = phy$edge.length[which(phy$tip.label == "Heracleum_sphondylium")],
                where = which(phy$tip.label == "Heracleum_sphondylium"))

# Make Spartina alterniflora, Spartina maritima, and Spartina x townsendii as sister to S. anglica
phy <- bind.tip(phy, "Spartina_alterniflora", 
                edge.length = phy$edge.length[which(phy$tip.label == "Spartina_anglica")],
                where = which(phy$tip.label == "Spartina_anglica"))

phy <- bind.tip(phy, "Spartina_maritima", 
                edge.length = phy$edge.length[which(phy$tip.label == "Spartina_anglica")],
                where = which(phy$tip.label == "Spartina_anglica"))

# Make Tephroseris integrifolia sister to all Senecio sp.
n <- getMRCA(phy, phy$tip.label[grep("Senecio", phy$tip.label)])

# Branch length is average of distance from Senecio MRCA to Senecio sp.
el1 <- sapply(nodes[grep("Senecio_", names(nodes))], nodeheight, tree = phy)
el2 <- nodeheight(tree = phy, node = n)

phy <- bind.tip(phy, "Tephroseris_integrifolia", edge.length = mean(el1 - el2), where = n)

# Make Viola hirta as sister to V. reichenbachiana (both in Viola sect. Viola)
phy <- bind.tip(phy, "Viola_hirta", 
                edge.length = phy$edge.length[which(phy$tip.label == "Viola_reichenbachiana")],
                where = which(phy$tip.label == "Viola_reichenbachiana"))

# Write second modified phylogeny
write.nexus(phy, file = paste0(path_proc_data, "/Lim_etal_2014_mod2.nex"))

# Check for species missing from modified phylogeny
notInPhy <- which(!(stomata$species %in% phy$tip.label |
                      stomata$acceptedname %in% phy$tip.label))
stomata[notInPhy, ]; length(notInPhy)

rm("edge.lengths", "el1", "el2", "l1", "l2", "n", "nodes", "notInPhy", "tips", "tmp")

# Run PATHd8 (Britton et al. 2002)
cat(sprintf("Sequence length = 4692;
            
            %s
            
            mrca: Azolla_filiculoides, Solanum_nigrum, maxage=454;
            mrca: Osmunda_regalis, Asplenium_onopteris, minage=354;
            mrca: Pinus_sylvestris, Solanum_nigrum, minage=310;
            mrca: Nymphaea_alba, Solanum_nigrum, minage=140;
            mrca: Nymphaea_alba, Solanum_nigrum, maxage=217;
            mrca: Nymphaea_alba, Nuphar_lutea, minage=115;
            mrca: Chelidonium_majus, Solanum_nigrum, fixage=121;
            mrca: Quercus_robur, Bryonia_dioica, minage=84;
            
            name of mrca: Azolla_filiculoides, Solanum_nigrum, name=tracheophytes;
            name of mrca: Osmunda_regalis, Asplenium_onopteris, name=monilophytes;
            name of mrca: Pinus_sylvestris, Solanum_nigrum, name=seedplants;
            name of mrca: Nymphaea_alba, Solanum_nigrum, name=angiosperms;
            name of mrca: Nymphaea_alba, Nuphar_lutea, name=nymphaceae;
            name of mrca: Chelidonium_majus, Solanum_nigrum, name=eudicots;
            name of mrca: Quercus_robur, Bryonia_dioica, name=curcibitales_fagales;
            
            ", write.tree(phy)), file = "PATHd8/pathd8_in")

# Run PATHd8
pathd8 <- str_c(getwd(), "/PATHd8") %>% str_replace_all(" ", "\\\\ ")
system(sprintf("cd %s; ./PATHd8 pathd8_in pathd8_out", pathd8))

# Import ultrametric tree
pathd8_out <- read_lines("PATHd8/pathd8_out")
phy <- read.tree(text = substr(pathd8_out[11], 14, nchar(pathd8_out[11])))

rm("pathd8_out")

# Write ultrametric phylogeny
write.nexus(phy, file = str_c(path_proc_data, "/Lim_etal_2014_final.nex"))

# Export objects to ms
export2ms(c("nNotInPhy1", "nNotInPhy2"))
