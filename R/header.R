rm(list = ls())
graphics.off()

compilte_pathd8 <- FALSE # changed to TRUE if you need to compule PATHd8
if (compilte_pathd8) {
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
library(vioplot)

# Directories
pathRawData <- "RawData"
pathProcData <- "ProcData"
pathR <- "R"
pathMS <- "ms"
pathFigures <- "ms/Figures"



