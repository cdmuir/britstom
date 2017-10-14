source("R/header.R")

# British Ecological Flora
stomata <- read_csv(str_c(path_raw_data, "/stomata.csv"), 
                    col_types = cols(
                      species = col_character(),
                      ab_density = col_double(),
                      ad_density = col_double()
                    ))

# Salisbury 1927
salisbury <- read_csv(str_c(path_raw_data, "/salisbury_1927.csv"),
                      col_types = cols(
                        species = col_character(),
                        ad_density = col_double(),
                        ab_density = col_double()
                      ))

# Taxize Salisbury dataset - this is done as a loop to prevent timing out
# tnr <- tnrs(query = salisbury$species[1], source = "iPlant_TNRS")
# for (i in 2:nrow(salisbury)) {
#   tnr %<>% bind_rows(tnrs(query = salisbury$species[i], source = "iPlant_TNRS"))
# }
# write_csv(tnr, path = str_c(path_raw_data, "/taxized_salisbury_1927.csv"))
tnr <- read_csv(str_c(path_raw_data, "/taxized_salisbury_1927.csv"),
                col_types = cols(
                  submittedname = col_character(),
                  acceptedname = col_character(),
                  sourceid = col_character(),
                  score = col_double(),
                  matchedname = col_character(),
                  authority = col_character(),
                  uri = col_character()
                ))

# When TNRS does not give an accepted name, I am using these based on searching The Plant List
# http://www.theplantlist.org/
# Plant List corrections
old_name <- c("Betula alba",            "Betula pubescens",         "Fagus sylvatica",
              "Ulmus montana",          "Corylus avellana",         "Crataegus oxyacanthoides",
              "Rubus corylifolius",     "Cephalanthera pallens",    "Conopodium denudatum",
              "Helleborus foetidus",    "Primula acaulis",          "Byronia dioica",
              "Geranium phaeum",        "Geranium striatum",        "Hieracium sylvaticum",
              "Hypericum pulchrum",     "Lathyrus macrorrhizus",    "Luzula forsteri",
              "Luzula pilosa",          "Melittis melissophyllum",  "Sison amomum",
              "Tilia parvifolia",       "Tilia vulgaris",           "Vicia sylvatica",
              "Hordeum sylvaticum")
new_name <- c("Betula pubescens",       "Betula pubescens",         "Fagus sylvatica",
              "Ulmus glabra",           "Corylus avellana",         "Crataegus laevigata",
              "Rubus plicatus",         "Cephalanthera longifolia", "Conopodium majus",
              "Helleborus foetidus",    "Primula vulgaris",         "Byronia dioica",
              "Geranium phaeum",        "Geranium versicolor",      "Hieracium murorum",
              "Hypericum pulchrum",     "Lathyrus linifolius",      "Luzula forsteri",
              "Luzula pilosa",          "Melittis melissophyllum",  "Sison amomum",
              "Tilia cordata",          "Tilia x europaea",         "Vicia sylvatica",
              "Hordelymus europaeus")

tnr$acceptedname[match(old_name, tnr$submittedname)] <- new_name
tnr %<>% mutate(species = submittedname)

# Changed one taxonomic name (does not affect trait or phylogeny data)
tnr$acceptedname %<>% str_replace("Arrhenatherum elatius var. elatius", 
                                  "Arrhenatherum elatius")

# tnr$acceptedname %>% str_split(" ") %>% sapply(length) %>% is_greater_than(1) %>% all()
salisbury %<>% join(tnr, by = "species")
rm(new_name, old_name, tnr)


# Find which Salisbury data are already in BEF data
already_in <- which(salisbury$submittedname %in% stomata$species |
                      salisbury$acceptedname %in% stomata$species)
x1 <- salisbury$submittedname %in% stomata$species
salisbury$species[x1] <- salisbury$submittedname[x1]
x2 <- (salisbury$acceptedname %in% stomata$species) & !x1
salisbury$species[x2] <- salisbury$acceptedname[x2]

tmp <- inner_join(stomata, salisbury, by = "species")

# Overlapping data from BEF and Salisbury are clearly taken from Salisbury, but with some errors. Hence, when they do not match, I have use the Salisbury data I entered myself.
with(tmp, plot(ab_density.x, ab_density.y))
with(tmp, plot(ad_density.x, ad_density.y))

x <- tmp %>%
  use_series(ab_density.x) %>%
  equals(tmp %>% use_series(ab_density.y)) %>%
  not() %>%
  which()

tmp$ab_density.x[x] <- tmp$ab_density.y[x]

x <- tmp %>%
  use_series(ad_density.x) %>%
  equals(tmp %>% use_series(ad_density.y)) %>%
  not() %>%
  which()

tmp$ad_density.x[x] <- tmp$ad_density.y[x]

tmp %<>% 
  dplyr::select(species, acceptedname, ab_density.x, ad_density.x) %>%
  set_colnames(c("species", "acceptedname", "ab_density", "ad_density"))

stomata %<>% 
  left_join(tmp, by = c("species", "ab_density", "ad_density"))

# Add data from Salisbury not in BEF
stomata %<>% 
  bind_rows(salisbury[!(x1 | x2), colnames(.)]) 

# Remove duplicates and stop if varietes are left
stomata %<>% filter(!duplicated(species))
stopifnot(str_detect(c(stomata$species, stomata$acceptedname), " var. ") %>%
            na.omit() %>% any() %>% not())

rm(already_in, salisbury, tmp, x, x1, x2)

# Calculate stomatal ratio
stomata %<>% mutate(sr_propAd = ad_density / (ab_density + ad_density),
                    sr_even = calc_sr_even(ab_density, ad_density))
stomata %<>% filter(!is.na(sr_propAd))

# Photosynthetic pathway
photo <- read_csv(str_c(path_raw_data, "/photo.csv"),
                  col_types = cols(
                    species = col_character(),
                    photo = col_character()
                  ))

# All missing taxa are C3 
# Not C4 based on Sage et al 2011 (JXB 62:9 3155-3169)
# Not CAM based on photosynthetic pathway of congeneric species in the data
missing_photo <- 
  c("Alchemilla vulgaris",  "Asarum europaeum", "Byronia dioica",
    "Hieracium sylvaticum", "Luzula maxima",    "Prunus cerasus",
    "Pteridium aquilinum",  "Pyrus communis",   "Pyrus torminalis",
    "Quercus sessiliflora", "Rhamnus frangula", "Ribes grossularia",
    "Rubus corylifolius",   "Tilia vulgaris")
photo %<>% bind_rows(data.frame(species = missing_photo, photo = "c3"))
stopifnot(all(stomata$species %in% photo$species |
                stomata$acceptedname %in% photo$species))
rm(missing_photo)

# PLANTATT (Hill et al 2004) for lifeform and Ellenberg light indicator values
plantatt <- read_csv(str_c(path_raw_data, "/plantatt.csv"),
                     col_types = cols(
                       .default = col_character(),
                       Chg = col_double(),
                       Hght = col_integer(),
                       Len = col_double(),
                       NBI = col_integer(),
                       SBI = col_integer(),
                       GB = col_integer(),
                       IR = col_integer(),
                       CI = col_integer(),
                       Tjan = col_double(),
                       Tjul = col_double(),
                       Prec = col_integer(),
                       L = col_integer(),
                       F = col_integer(),
                       R = col_integer(),
                       N = col_integer(),
                       S = col_integer()
                     )) %>%
  dplyr::select(species = `Taxon name`, 
                lifeform = LF1, perennation1 = P1, perennation2 = P2, woodiness = W, 
                height = Hght, ellenberg_light = L)

# Define habit
plantatt$habit <- factor(NA, levels = c("annual", "biennial", "perennial",
                                        "shrub", "tree"))
  
# Annual if P1 == annual and is.na(P2)
x <- which(plantatt$perennation1 == "a" & is.na(plantatt$perennation2))
plantatt$habit[x] <- "annual"

# Biennial if P1 == biennial & is.na(P2)
# or P1 == annual & !is.na(P2)l
# or P1 == biennial & !is.na(P2)
# or P1 == perennial & !is.na(P2)

x <- which(plantatt$perennation1 == "b" & is.na(plantatt$perennation2) |
             plantatt$perennation1 == "a" & !is.na(plantatt$perennation2) |
             plantatt$perennation1 == "b" & !is.na(plantatt$perennation2) |
             plantatt$perennation1 == "p" & !is.na(plantatt$perennation2))

plantatt$habit[x] <- "biennial"

# Perennial if P1 == perennial and is.na(P2)
x <- which(plantatt$perennation1 == "p" & is.na(plantatt$perennation2))
plantatt$habit[x] <- "perennial"

# Reclassify as shrub if woody
x <- which(plantatt$woodiness == "w")
plantatt$habit[x] <- "shrub"

# Reclassify as tree if height > 4 m
shrub_cutoff <- 4
x <- which(plantatt$woodiness == "w" & plantatt$height >= (shrub_cutoff * 1e2))
plantatt$habit[x] <- "tree"

rm(x)
plantatt %<>% dplyr::select(species, lifeform, ellenberg_light, habit)

# stomata$species[!(stomata$species %in% plantatt$species |
#                    stomata$acceptedname %in% plantatt$species)] # Missing data for several

# Combine bulbous with nonbulbous geophytes, hydrophytes with annual hydrophytes,
# and phanerophytes with nanophanerophytes
plantatt$lifeform %<>% mapvalues(from = c("Gb", "Hz", "Pn"), 
                                 to = c("Gn", "Hy", "Ph")) 

# Join datasets
xp <- match(stomata$species, photo$species)
xpa <- match(stomata$acceptedname, photo$species)
xp[is.na(xp)] <- xpa[is.na(xp)]

xa <- match(stomata$species, plantatt$species)
xaa <- match(stomata$acceptedname, plantatt$species)
xa[is.na(xa)] <- xaa[is.na(xa)]

stomata$photo <- photo$photo[xp]
stomata$lifeform <- plantatt$lifeform[xa]
stomata$ellenberg_light <- plantatt$ellenberg_light[xa]
stomata$habit <- plantatt$habit[xa]

# Remove missing values
stomata %<>% filter(!is.na(sr_propAd), !is.na(lifeform), !is.na(ellenberg_light))

# Removing species that are classified as both c3 and cam
stomata %<>% filter(!(species %in% c("Sedum acre", "Sedum telephium")))
stopifnot(photo$species %>% .[duplicated(.)] %>% is_in(stomata$species) %>% any() %>% not())
stopifnot(photo$species %>% .[duplicated(.)] %>% is_in(stomata$acceptedname) %>% any() %>% not())

rm(photo, plantatt, xa, xaa, xp, xpa)

# Export
stomata$species %<>% str_replace_all(" ", "_")
stomata$acceptedname %<>% str_replace_all(" ", "_")
write_csv(stomata, str_c(path_proc_data, "/stomata.csv"))

# Export objects to ms
stomata_unfiltered <- stomata
export2ms(c("shrub_cutoff", "stomata_unfiltered"))
