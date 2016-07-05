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




### pollen transport analysis

# summarise pollen transport for each insect species within each sample
# this aggregates all pollen grains carried by an insect species for each plant species within each sample

dframe2 <- ddply(dframe1r, .(Family_Species,Sample), numcolwise(sum))

# create a suitable matrix from the dframe
dframe2m <- ddply(dframe2, .(Family_Species), numcolwise(sum))
rownames(dframe2m) <- dframe2m[,1]
dframe2m <- dframe2m[,-1]

### plot a full-system network
plotweb(dframe2m)
# networklevel(dframe3m)
# trying to run networklevel() on this full network causes my computer to crash, so I've concealed it behind a # to stop it from running every time!


### that's fine, but we want networks (and descriptors) for each sample

# split each sample into a separate dataframe
dframesPT <- split(dframe2, list(dframe2$Sample))  # this creates a list of smaller dframes, one for each level of sample
summary(dframesPT)
# for example...
summary(dframesPT[1]) # the first dframe in the list is site F1 on sampling day 10
dframesPT[1]

### using a loop and a custom function called 'network', ...

### prepare the dataframes for network analysis and run it (this command will take a few minutes to run)
# the command will print "Fail" each time it detects a dataframe with too little data to produce a network
metricsPT <- lapply(dframesPT, network)
summary(metricsPT)

###

### we now have a list of dataframes, each one containing the network metricsPT of one site

# for sites that had too little data, the dataframes just contain the value "Fail" as a character
# combine all dataframes together
metricsPT.merge <- do.call("cbind", metricsPT)    # merge the data with one sample per column
colnames(metricsPT.merge) <- names(metricsPT)     # assign the sample names to each column
metricsPT.merge <- data.frame(t(metricsPT.merge))
metricsPT.merge

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
metricsPT.full <- merge(dframeP, metricsPT.merge, by=0)


### before we can analyse this data, we need to remove any samples which have failed

# subset the columns to keep only those which do not have "Fail" values
metricsPT.good <- metricsPT.full[ which( ! metricsPT.full$connectance %in% "Fail") , ]

# remove the duplicate 'row.names' column - this information is the same as 'Sample'
metricsPT.good <- metricsPT.good[,-1]


### hooray! we now have network metrics for each sample with sufficient data to produce a network


### output the dataframe to a .txt file to use in downstream analysis
write.table(metricsPT.good, "Data\\PTNetworkmetrics.txt", sep="\t", row.names=FALSE)










### flower-visitor analysis

# change each interaction to a 1 for flower-visitor analysis
dframe1f <- dframe1r
dframe1f[,3:73][dframe1f[,3:73] > 0] <- 1


# summarise plant-insect interactions for each insect species within each sample
# this produces semi-quantitative data - e.g. if two individuals of an insect species
# interact with a plant species in a sample,
# the entry for that interaction within that sample will be 2

dframe3 <- ddply(dframe1f, .(Family_Species,Sample), numcolwise(sum))
  
# create a suitable matrix from the dframe
dframe3m <- ddply(dframe3, .(Family_Species), numcolwise(sum))
rownames(dframe3m) <- dframe3m[,1]
dframe3m <- dframe3m[,-1]

### plot a full-system network
plotweb(dframe3m)
# networklevel(dframe3m)
# trying to run networklevel() on this full network causes my computer to crash, so I've concealed it behind a # to stop it from running every time!


### that's fine, but we want networks (and descriptors) for each sample

# split each sample into a separate dataframe
dframesFV <- split(dframe3, list(dframe3$Sample))  # this creates a list of smaller dframes, one for each level of sample
summary(dframesFV)
# for example...
summary(dframesFV[1]) # the first dframe in the list is site F1 on sampling day 10
dframesFV[1]

### using a loop and a custom function called 'network', ...

### prepare the dataframes for network analysis and run it (this command will take a few minutes to run)
# the command will print "Fail" each time it detects a dataframe with too little data to produce a network
metricsFV <- lapply(dframesFV, network)
summary(metricsFV)

###

### we now have a list of dataframes, each one containing the network metricsFV of one site

# for sites that had too little data, the dataframes just contain the value "Fail" as a character
# combine all dataframes together
metricsFV.merge <- do.call("cbind", metricsFV)    # merge the data with one sample per column
colnames(metricsFV.merge) <- names(metricsFV)     # assign the sample names to each column
metricsFV.merge <- data.frame(t(metricsFV.merge))
metricsFV.merge

### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the previously used sample-level pollen dataframe:
# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
metricsFV.full <- merge(dframeP, metricsFV.merge, by=0)


### before we can analyse this data, we need to remove any samples which have failed

# subset the columns to keep only those which do not have "Fail" values
metricsFV.good <- metricsFV.full[ which( ! metricsFV.full$connectance %in% "Fail") , ]
  
# remove the duplicate 'row.names' column - this information is the same as 'Sample'
metricsFV.good <- metricsFV.good[,-1]



### output the dataframe to a .txt file to use in downstream analysis
write.table(metricsFV.good, "Data\\FVNetworkmetrics.txt", sep="\t", row.names=FALSE)




