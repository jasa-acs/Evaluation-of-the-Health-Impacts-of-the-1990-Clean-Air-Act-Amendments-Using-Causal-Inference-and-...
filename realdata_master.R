library(mvtnorm)
library(dplyr)
library(splines)
library(BayesTree)
library(xtable)
library(ggplot2)
library(mgcv)

## create dataset ##
source('create_dataset.R')

# clear the workspace using
#rm(list=ls())

## run the analysis ##
source('analysis_wholeUS.R')

# clear the workspace using
#rm(list=ls())

## make table 2, export to file 'table_2.txt' ##
source('table_2.R')

# clear the workspace using
#rm(list=ls())

## make figure 2, export to file 'figure_2.pdf' ##
source('figure_2.R')
