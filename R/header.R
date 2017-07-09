rm(list = ls())
graphics.off()

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

# Directories
pathRawData <- "RawData"
pathProcData <- "ProcData"
pathR <- "R"
pathMS <- "ms"
pathFigures <- "Figures"



