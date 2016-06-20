##################################################
####   Script for network analysis by sample  ####
##################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","car","ggplot2","RVAideMemoire","arm","bipartite","plyr")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded


### load up Callum's custom set of functions
source("CheckResidsFunction.R")


### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/MatrixNoct.txt", header=TRUE)

summary(dframe1) # Check it's imported correctly

### prepare the data for network analysis

# trim off extra columns
dframe1r <- dframe1[,c(2,9:82)]

# change each interaction to a 1
dframe1r[,3:75][dframe1r[,3:75] > 0] <- 1


# summarise plant-insect interactions for each insect species within each sample
# this produces semi-quantitative data - e.g. if two individuals of an insect species
# interact with a plant species in a sample,
# the entry for that interaction within that sample will be 2

dframe2 <- ddply(dframe1r, .(Family_Species,Sample), numcolwise(sum))
  
# create a suitable matrix from the dframe
dframe2m <- ddply(dframe2, .(Family_Species), numcolwise(sum))
rownames(dframe2m) <- dframe2m[,1]
dframe2m <- dframe2m[,-1]

### plot a full-system network
plotweb(dframe2m)
networklevel(dframe2m)


### that's fine, but we want networks (and descriptors) for each sample

# split each sample into a separate dataframe
dframes <- split(dframe2, list(dframe2$Sample))

summary(dframes)

# using a loop, ...

# summarise each dataframe
lapply(dframes, function(x) summary(x))

# remove the Sample column
dframes1 <- lapply(dframes, function(x) x[,c(1,3:75)])


plotweb(dframes1[1])





# make the first column the row names
dframes2 <- lapply(dframes1, function(x) 
                      rownames(x) <- x[,1])

dframes3 <- lapply(dframes2, function(x)
                x <- x[,-1])
  
  
  samp2 <- samp[,-1]
       rownames(samp2) <- samp[,1])

lapply(dframes, function(x)
  plotweb(x)
)


####





plotweb(CALLUM)
#change arrows
plotweb(CALLUM, arrow="down")
#change colours etc
plotweb(CALLUM, arrow="down", col.high = "red")
#change individual colours
plotweb(CALLUM, arrow="down", col.high = c("red", "green", "yellow", "red", "blue", "black", "blue", "red"))
# colours above and below
plotweb(CALLUM, arrow="down", col.high = c("red", "green", "yellow", "red", "blue", "black", "blue", "red"), col.low = c("green","yellow", "red", "blue", "black", "blue", "red", "brown"))
# see col.interaction for changing link colours and use c() as above for individuals
#calculate network metrics
networklevel(CALLUM)
