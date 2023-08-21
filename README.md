# ViralBottleneck
## Description
This package is used for estimating viral transmission bottleneck size using different methods.
The main steps to use program are:
1. Transmission object creation
2. Transmission bottleneck size calculation
   
There are the he [manual](manual_and_tutorial/ViralBottleneck_manual_0.1.0.pdf) and [turtorial](manual_and_tutorial/Tutorial.pdf).

The contents in tables in [turtorial](manual_and_turtorial/Tutorial.pdf) are little crowded. Therefore, there are also a html file of turtorial in [manual_and_turtorial folder](manual_and_tutorial). It would be better for users to read. Please download html file and open in the browser. The test data for "ViralBottleneck" is in [test_dataset folder](test_dataset).

## Download
The codes for downloading package "ViralBottleneck" are below: 
```
library(devtools)
install_github("BowenArchaman/ViralBottleneck")
library(ViralBottleneck)
```

## Requirements
- R 4.2.2
- methods
- utils
- ggplot2
- pbapply
- rmutil
