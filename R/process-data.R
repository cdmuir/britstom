source("R/header.R")

# Stomata

stomata <- read_csv(str_c(pathRawData, "/stomata.csv"))
stomata %<>% mutate(sr_propAd = ad_density / (ab_density + ad_density),
                    sr_even = min(c(ab_density, ad_density)) / 
                      max(c(ab_density, ad_density)))
stomata %<>% filter(!is.na(sr_propAd))

# Photosynthetic pathway
photo <- read_csv(str_c(pathRawData, "/photo.csv"))
all(stomata$species %in% photo$species) # should be true

# PLANTATT (Hill et al 2004) for lifeform and Ellenberg light indicator values
plantatt <- read_csv(str_c(pathRawData, "/plantatt.csv")) %>%
  select(species = `Taxon name`, lifeform = LF1, ellenberg_light = L, 
         height = Hght, wh = W)
# stomata$species[which(!stomata$species %in% plantatt$species)] # Missing data for a few

# Combine bulbous with nonbulbous geophytes, hydrophytes with annual hydrophytes,
# and phanerophytes with nanophanerophytes
plantatt$lifeform %<>% mapvalues(from = c("Gb", "Hz", "Pn"), to = c("Gn", "Hy", "Ph")) 

# Join datasets
stomata %<>% join(photo, by = "species")
stomata %<>% join(plantatt, by = "species")

# Remove missing values
stomata %<>% filter(!is.na(sr_propAd), !is.na(lifeform), !is.na(ellenberg_light))

# Removing species that are classified as both c3 and cam
stomata %<>% filter(!(species %in% c("Sedum acre", "Sedum telephium")))

# Export  
write_csv(stomata, str_c(pathProcData, "/stomata.csv"))

# OLD - rm("photo", "plantatt", "xp", "xa")

# OLD - Ellenberg light index
# ellenberg_light <- read.csv(paste0(pathRawData, "/ellenberg_light.csv"), stringsAsFactors = F, 
#                            strip.white = T)

# OLD - Life-form
# lifeform <- read.csv(paste0(pathRawData, "/life-form.csv"), stringsAsFactors = F, strip.white = T)

# x1 <- match(stomata$species, ellenberg_light$species)
# x2 <- match(stomata$species, lifeform$species)
# stomata$ellenberg_light <- ellenberg_light$ellenberg_light[x1]
# stomata$lifeform <- lifeform$lifeform[x2]
