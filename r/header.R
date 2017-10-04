rm(list = ls())
graphics.off()

# remove before final publication
# file.copy("~/Google Drive/CommonResources/refs.bib",
#           "~/Google Drive/britstom/ms/refs.bib", overwrite = TRUE)

compile_pathd8 <- FALSE # changed to TRUE if you need to compile PATHd8
if (compile_pathd8) {
  system("cd PATHd8; cc PATHd8.c -O3 -lm -o PATHd8")
}

# Libraries
library(ape)
library(dplyr)
library(lavaan)
library(magrittr)
library(phylolm)
library(phytools)
library(plyr)
library(readr)
library(rncl)
library(RColorBrewer)
library(Rphylopars)
library(stringr)
library(taxize)
library(tidyr)
library(vioplot)

# Directories
path_raw_data <- "raw-data" %>% normalizePath()
path_proc_data <- "proc-data" %>% normalizePath()
path_r <- "r" %>% normalizePath()
path_objects <- str_c(path_r, "/objects") %>% normalizePath()
path_ms <- "ms" %>% normalizePath()
path_figures <- str_c(path_ms, "/figures") %>% normalizePath()

source(str_c(path_r, "/functions.R"))
