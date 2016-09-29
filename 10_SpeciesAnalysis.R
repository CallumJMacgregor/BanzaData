########################################################
####   Script for species-level analysis by season  ####
########################################################


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

dframe1$SeasonTreatment <- do.call(paste, c(dframe1[c("Season","Year","Treatment")], sep = "_"))  

# select columns in specific order: Species, Network ID, Pollen
dframe1r <- dframe1[,c(2,9,13:80)]
dframe1s <- dframe1[,c(2,81,13:80)]



### pollen transport analysis

### per sample networks

# summarise pollen transport for each insect species within each sample
# this aggregates all pollen grains carried by an insect species for each plant species within each sample

dframe2 <- ddply(dframe1r, .(Family_Species,Sample), numcolwise(sum))

# create a suitable matrix from the dframe
dframe2m <- ddply(dframe2, .(Family_Species), numcolwise(sum))
rownames(dframe2m) <- dframe2m[,1]
dframe2m <- dframe2m[,-1]

### plot a full-system network
plotweb(dframe2m)
# trying to run networklevel() on this full network causes my computer to crash, so I've concealed it behind a # to stop it from running every time!


### that's fine, but we want networks (and descriptors) for each season
### per season networks

# summarise pollen transport for each insect species within each season
# this aggregates all pollen grains carried by an insect species for each plant species within each sample

dframe3 <- ddply(dframe1s, .(Family_Species,SeasonTreatment), numcolwise(sum))

### that's fine, but we want networks (and descriptors) for each sample

# split each sample into a separate dataframe
dframesPTS <- split(dframe3, list(dframe3$SeasonTreatment))  # this creates a list of smaller dframes, one for each level of sample
summary(dframesPTS)
# for example...
summary(dframesPTS[1]) # the first dframe in the list is site F1 on sampling day 10
dframesPTS[1]


### using a loop and a custom function called 'specieslev', ...

### prepare the dataframes for species-level network analysis and run it (this command will take a few minutes to run)
# the command will print "Fail" each time it detects a dataframe with too little data to produce a network
metricsPTS <- lapply(dframesPTS, specieslev)
summary(metricsPTS)


###

### we now have a list of dataframes, each one containing the network metricsPT of one site

# for sites that had too little data, the dataframes just contain the value "Fail" as a character
# combine all dataframes together
metricsPTS.merge <- do.call("rbind", metricsPTS)    # tack the data together in a single dframe


### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the previously used sample-level pollen dataframe:
dframeP<-read.table("Data/SamplesNoct.txt", header=TRUE)
summary(dframeP)

# we'll need a bit more data later on, so save this in a separate dframe
dframePS <- dframeP[,1:8]

# we'll need to replicate the SeasonTreatment variable
dframePS$sample <- do.call(paste, c(dframePS[c("Season","Year","Treatment")], sep = "_")) 
dframePS$sample <- factor(dframePS$sample)

# we need to collapse this down into the same 18 rows as the metrics
dframePSc <- ddply(dframePS, .(Treatment,Season,Year,sample), numcolwise(sum))

# we don't need the SamplingDay variable any more
dframePSc <- dframePSc[,-5]

# now merge them together (there is only one matching colname, "sample", which is the one we want to merge by...
# ...so we don't need to specify)
metricsPTS.full <- merge(dframePSc,metricsPTS.merge)

metricsPTS.full <- metricsPTS.full[,c(length(metricsPTS.full),(length(metricsPTS.full)-1),1:(length(metricsPTS.full)-2))]


### we want to analyse this data separately for insects and plants

# subset the columns to keep only those which do not have "Fail" values
metricsPTS.plants <- metricsPTS.full[ which( ! metricsPTS.full$level %in% "insect") , ]
metricsPTS.insects <- metricsPTS.full[ which( ! metricsPTS.full$level %in% "plant") , ]


### hooray! we now have species-level metrics for each season


### output the dataframe to a .txt file to use in downstream analysis
write.table(metricsPTS.full, "Data\\SpecieslevelmetricsPTbySeason.txt", sep="\t", row.names=FALSE)




