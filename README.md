# britstom
Light and growth form interact to shape stomatal ratio among British angiosperms

This project was developed by [Chris Muir](www.chrisdmuir.com).

More information about the method is available in a preprint which you can find on [biorxiv](http://biorxiv.org/?????) or on [github](https://github.com/cdmuir/ms/ms.pdf).

## Downloading data and code 

1. Download or clone this repository to your machine.
2. Open `britstom.Rproj` in [RStudio](https://www.rstudio.com/)

## Generating manuscript

You can source all the code you need in the correct order using `r/run-all.R`. Even if you don't want to run all the code, you may need to install some packages (`r/install-packages.R`) and attach them (`r/header.R`).

- To use premade R output, simply open `ms/ms.Rnw` and compile using RStudio.
- To rerun all analyses, you will first need to download and compile [PATHd8](http://www2.math.su.se/PATHd8/) locally in the PATHd8 directory associated with this project. To compile, change code in `r/header.R` to:

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
