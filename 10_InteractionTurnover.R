##################################################
####   Script for network analysis by sample  ####
##################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("betalink","plyr")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded


### load up Callum's custom set of functions
f <- c("NetworkFunction.R")
lapply(f, source)



### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/MatrixNoct.txt", header=TRUE)

summary(dframe1) # Check it's imported correctly

### prepare the data for network analysis

# add a column to identify which network each row belongs to

dframe1$Network <- do.call(paste, c(dframe1[c("Season","Year","Treatment")], sep = "_"))
dframe1$Network <- as.factor(dframe1$Network)
summary(dframe1$Network)

# and another for which network pair each row belongs to

dframe1$Pair <- do.call(paste, c(dframe1[c("Season","Year")], sep = "_"))
dframe1$Pair <- as.factor(dframe1$Pair)
summary(dframe1$Pair)

# trim off extra columns
dframe1r <- dframe1[,c(2,7,(length(dframe1)-1),length(dframe1),13:(length(dframe1)-2))]

# change each interaction to a 1 for flower-visitation
dframe1r[,5:length(dframe1r)][dframe1r[,5:length(dframe1r)] > 0] <- 1

# summarise plant-insect interactions for each insect species within each sample
# this produces semi-quantitative data - e.g. if two individuals of an insect species
# interact with a plant species in a sample,
# the entry for that interaction within that sample will be 2

dframe2 <- ddply(dframe1r, .(Family_Species,Treatment,Pair,Network), numcolwise(sum))

# create separate dframes depending on treatment

dframe.fire <- dframe2[ which(dframe2$Treatment=="Fire"), ]
dframe.nofire <- dframe2[ which(dframe2$Treatment=="NoFire"), ]


# split each sample into a separate dataframe
fires <- split(dframe.fire, list(dframe.fire$Network))  # this creates a list of smaller dframes, one for each level of sample
summary(fires)

nofires <- split(dframe.nofire, list(dframe.nofire$Network))
summary(nofires)



network.pairs <- mapply(c, fires, nofires, SIMPLIFY = FALSE)
summary(network.pairs)


# this is not yet in proper matrix form for network analysis, so use custom function 'prepare':
dframes.prep <- lapply(dframes, prepare)









########################### development #######################


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





