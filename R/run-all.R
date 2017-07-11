# source("R/install-packages.R")
# source("R/update-packages.R")

# source("R/session-info.R")

# Prepare data
source("R/01-process-data.R")
source("R/02-make-phylogeny.R")
source("R/03-plot-lifeform.R")
source("R/04-filter-data.R") # starting here (if not before) need place to export stuff for use in ms
source("R/05-modify-phylogeny.R")

# Analyze data
source("R/06-analyse-light-vs-lifeform.R")
source("R/07-analyse-stomata.R")

# NEED TO GET PATH ANALYSIS STUFF
# NEED TO KNIT REPORT