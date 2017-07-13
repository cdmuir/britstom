# britstom
Light and life form interact to shape stomatal ratio among British angiosperms

Biologists are increasingly using curated, public data sets to conduct phylogenetic comparative analyses. Unfortunately, there is often a mismatch between species for which there is phylogenetic data and those for which other data is available. As a result, researchers are commonly forced to either drop species from analyses entirely or else impute the missing data.

In this package we have implemented a simple solution to increase the overlap while avoiding potential the biases introduced by imputing data.  If some external topological or taxonomic information is available, this can be used to maximize the overlap between the data and the phylogeny. The algorithms in `phyndr` replace a species lacking data with a species  that has data. This swap can be made because for those two species, all phylogenetic relationships are exactly equivalent.

This project was developed by [Chris Muir](www.chrisdmuir.com).

More information about the method is available in a preprint which you can find on [biorxiv](http://biorxiv.org/?????) or on [github](https://github.com/traitecoevo/ms/ms.pdf).

## Downloading data and code 

1. Download or clone this repository to your machine.
2. Open `britstom.Rproj` in [RStudio](https://www.rstudio.com/)

## Generating manuscript

You can source all the code you need in the correct order using `r/run-all.R`. Even if you don't want to run all the code, you may need to install some packages (`r/install-packages.R`) and attach them (`r/header.R`).

- To use premade R output, simply open `ms/ms.Rnw` and compile using RStudio.
- To rerun all analyses, first compile locally [PATHd8](http://www2.math.su.se/PATHd8/). Change code in `r/header.R` to:

```
compile_pathd8 <- TRUE
if (compile_pathd8) {
  system("cd PATHd8; cc PATHd8.c -O3 -lm -o PATHd8")
}

```
- Next, source `r/run-all.R` in the R Console:

```

# This will take several minutes to run
source("r/run-all.R")

```

[![CC BY](http://i.creativecommons.org/l/by/3.0/88x31.png)](http://creativecommons.org/licenses/by/3.0/)
