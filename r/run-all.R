# source("r/install-packages.R")

# source("r/session-info.R")

source("r/header.R")

# Prepare data
source("r/01-process-data.R")
source("r/02-make-phylogeny.R")
source("r/03-plot-lifeform.R")
source("r/04-filter-data.R")
source("r/05-modify-phylogeny.R")

# Analyze data
source("r/06-analyse-light-vs-lifeform.R")
source("r/07-analyse-stomatal-ratio.R")
source("r/08-analyse-stomatal-density.R")
source("r/09-analyse-phyloSEM.R")
