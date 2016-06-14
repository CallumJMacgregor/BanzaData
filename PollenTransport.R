########################################################
####   Script for basic pollen transport analysis   ####
########################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","car","ggplot2")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


library(lme4)
library(car)
library(ggplot2)

### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/MatrixNoct.txt", header=TRUE)

summary(dframe1) # Check it's imported correctly

### tell R that SampleID and SlideNumber should be treated as factors
dframe1$SampleID <- factor(dframe1$SampleID)
dframe1$SlideNumber <- factor(dframe1$SlideNumber)

summary(dframe1)

### total up the pollen grains for each insect
dframe1$PollenLoad<-rowSums(dframe1[c(7:79)])

### total up the number of pollen species for each insect (i.e. how many columns do not contain 0)
dframe1$PollenTypes<-rowSums(dframe1[c(7:79)] != 0)

### now you're ready to start looking for patterns!