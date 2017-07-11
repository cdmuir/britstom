rm(list = ls())
graphics.off()

compile_pathd8 <- FALSE # changed to TRUE if you need to compile PATHd8
if (compile_pathd8) {
  system("cd PATHd8; cc PATHd8.c -O3 -lm -o PATHd8")
}

source("R/functions.R")

# Libraries
library(ape)
library(dplyr)
library(magrittr)
library(phylolm)
library(phytools)
library(plyr)
library(readr)
library(rncl)
library(Rphylopars)
library(rstanarm)
library(stringr)
library(taxize)
library(vioplot)

# Directories
path_raw_data <- "raw-data"
path_proc_data <- "proc-data"
path_r <- "r"
  path_objects <- "r/objects"
pathM_ms <- "ms"
  path_figures <- "ms/figures"



