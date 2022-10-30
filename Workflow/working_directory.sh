#!/bin/bash

#The path / working directory you choose
userInput=$1

#Select which scripts to run, add "#" before each line if want to skip it.
#But note that these scripts should be run based on their order (except that 2, 3 and 4 can be run at the same time).

Rscript ./1_Sobol.R "$userInput"

Rscript ./2_Kriging.R "$userInput"

Rscript ./3_BASS.R "$userInput"

Rscript ./4_AKMCS.R "$userInput"

Rscript ./5_traceplot_data.R "$userInput"

Rscript ./6_traceplot.R "$userInput"

Rscript ./7_Gridplot.R "$userInput"

Rscript ./8_Supplementary.R "$userInput"
