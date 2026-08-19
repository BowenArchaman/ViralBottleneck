# ViralBottleneck
## Description
This package is used for estimating viral transmission bottleneck size using different methods.
The main steps to use the program are:
1. Create the transmission object
2. Calculate the transmission bottleneck size 
   
There are the links for [manual](manual_and_tutorial/ViralBottleneck_0.1.0.pdf) and [tutorial](manual_and_tutorial/ViralBottleneck_Tutorial.pdf).

The test data for "ViralBottleneck" is in [test_dataset folder](test_dataset).

## Download
To install the package in R use the following commands: 
```
library(devtools)
install_github("BowenArchaman/ViralBottleneck", build_vignettes = TRUE)
library(ViralBottleneck)
```

## Requirements
- R 4.2.2
- methods
- utils
- ggplot2
- pbapply
- rmutil

## Citation
Please cite the following manuscript when you use the R package:

Zheng B, Johnson PCD, Hughes J. **ViralBottleneck: an R package for estimating viral transmission bottlenecks from deep sequencing data using multiple methods.** Virus Evol. 2025 Sep 19;11(1):veaf071. [doi: 10.1093/ve/veaf071](https://doi.org/10.1093/ve/veaf071). Erratum in: Virus Evol. 2026 Aug 14;12(1):veag055. doi: 10.1093/ve/veag055. PMID: 41089325; PMCID: PMC12516950.
