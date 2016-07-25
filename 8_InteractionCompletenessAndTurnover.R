################################################################
####   Script for network analysis by season and treatment  ####
################################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("betalink","plyr")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded


### load up Callum's custom set of functions
f <- c("NetworkFunction.R","SpecAccumFunctions.R")
lapply(f, source)



### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/MatrixNoct.txt", header=TRUE)

summary(dframe1) # Check it's imported correctly

dframe1$ExactSeason <- do.call(paste, c(dframe1r[c("Season","Year","Treatment")], sep = "_"))
dframe1$ExactSeason <- as.factor(dframe1$ExactSeason)


### prepare the data for network analysis

# trim off extra columns
dframe1r <- dframe1[,c(2,7,10,12:length(dframe1))]

# change each interaction to a 1
dframe1r[,5:length(dframe1r)][dframe1r[,3:length(dframe1r)] > 0] <- 1



# summarise plant-insect interactions for each insect species within each Season-Year-Treatment
# this produces semi-quantitative data - e.g. if two individuals of an insect species
# interact with a plant species in a sample,
# the entry for that interaction within that sample will be 2

dframe2 <- ddply(dframe1r, .(Family_Species,Season,Year,Treatment,ExactSeason), numcolwise(sum))


# split each sample into a separate dataframe
dframes <- split(dframe2, list(dframe2$ExactSeason))  # this creates a list of smaller dframes, one for each level of sample
summary(dframes)


# first let's use another of Callum's custom functions to plot all the species accumulation curves
lapply(dframes, sacplot, cols=5)















# create a vector for the number of samples
x <- c(1:length(dframes))
x

# create a list of possible pairs of dframes
pairs <- combn(x,2)







# this is not yet in proper matrix form for network analysis, so use custom function 'prepare':
dframes.prep <- lapply(dframes, prepare)

# and now use 'prepare_networks' from the betalink package to convert these into the correct format for betalink (igraph)
networks <- prepare_networks(dframes.prep)
summary(networks)

data.frame(B-div) <- network_betadiversity(networks)










graph1 <- networks[1]
plot(graph1)

# measure dissimilarity between networks; there are 23 options for B-diversity measures but let's use the default
distance <- beta_os_prime(networks)
distance

B-div <- network_betadiversity(networks)

plot <- network_betaplot(networks[1], networks[2], na="grey20", nb="black",ns="grey70")












