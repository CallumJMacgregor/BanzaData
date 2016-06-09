#########################################
####   Script for preparing data for analysis   ####
#########################################


### Clear the workspace
rm(list=ls())


### load the libraries you need
j <- c("reshape2","ggplot2")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


library(reshape2)
library(ggplot2)

