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
f <- c("NetworkFunction.R")
lapply(f, source)



### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/MatrixNoct.txt", header=TRUE)

summary(dframe1) # Check it's imported correctly

### prepare the data for network analysis

# trim off extra columns
dframe1r <- dframe1[,c(2,9:80)]

# change each interaction to a 1 for flower-visitor analysis
dframe1r[,3:73][dframe1r[,3:73] > 0] <- 1


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
# networklevel(dframe2m)
# trying to run networklevel() on this full network causes my computer to crash, so I've concealed it behind a # to stop it from running every time!


### that's fine, but we want networks (and descriptors) for each sample

# split each sample into a separate dataframe
dframes <- split(dframe2, list(dframe2$Sample))  # this creates a list of smaller dframes, one for each level of sample
summary(dframes)
# for example...
summary(dframes[1]) # the first dframe in the list is site F1 on sampling day 10
dframes[1]

### using a loop and a custom function called 'network', ...

### prepare the dataframes for network analysis and run it (this command will take a few minutes to run)
# the command will print "Fail" each time it detects a dataframe with too little data to produce a network
metrics <- lapply(dframes, network)
summary(metrics)

###

### we now have a list of dataframes, each one containing the network metrics of one site

# for sites that had too little data, the dataframes just contain the value "Fail" as a character
# combine all dataframes together
metrics.merge <- do.call("cbind", metrics)    # merge the data with one sample per column
colnames(metrics.merge) <- names(metrics)     # assign the sample names to each column
metrics.merge <- data.frame(t(metrics.merge))
metrics.merge

### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the previously used sample-level pollen dataframe:
dframeP<-read.table("Data/SamplesNoct.txt", header=TRUE)
summary(dframeP)

# we don't need any of the pollen data, so:
dframeP <- dframeP[,1:3]

# but we do need a Treatment variable, so:
dframeP$Treatment <- ifelse(grepl("NF",dframeP$Site),"NoFire","Fire")
dframeP$Treatment <- factor(dframeP$Treatment)

# for merge to work, we need to set the row names
rownames(dframeP) <- dframeP$Sample

# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
metrics.full <- merge(dframeP, metrics.merge, by=0)


### before we can analyse this data, we need to remove any samples which have failed

# subset the columns to keep only those which do not have "Fail" values
metrics.good <- metrics.full[ which( ! metrics.full$connectance %in% "Fail") , ]
  
# remove the duplicate 'row.names' column - this information is the same as 'Sample'
metrics.good <- metrics.good[,-1]


### hooray! we now have network metrics for each sample with sufficient data to produce a network


### output the dataframe to a .txt file to use in downstream analysis
write.table(metrics.good, "Data\\NetworkMetrics.txt", sep="\t", row.names=FALSE)




