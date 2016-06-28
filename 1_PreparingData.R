####################################################
####   Script for preparing data for analysis   ####
####################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("reshape2","plyr")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded


### read in the data - this is the raw data as input by Paula. 
# I have manually added a column called 'SampleID' with a unique number for each row of the dataframe;
# I have renamed "Family/Species" to "Family_Species", 
# and added _1 to the first PollenType and PollenNumber column for consistency
# Finally, I've put a 0 in "PollenNumber_1" for insects with no pollen - this is important
dframe1<-read.csv("Data/BanzaDataNoct.csv", header=TRUE)

summary(dframe1) # Check it's imported correctly

### Remove redundant PollenTypes 14 & 15
dframe1 <- dframe1[,1:32]
summary(dframe1)

### Melt the dataframe
names(dframe1) # check what columns you have and what they're called

# create a dframe with all the family/species names for each sample
melt1 <- melt(dframe1[,c("SampleID","Family_Species","Date","Site","SlideNumber","PollenCount",
                         "PollenType_1","PollenType_2","PollenType_3","PollenType_4",
                         "PollenType_5","PollenType_6","PollenType_7","PollenType_8",
                         "PollenType_9","PollenType_10","PollenType_11","PollenType_12","PollenType_13")],
                      id=c("SampleID","Family_Species","Date","Site","SlideNumber","PollenCount"))

#label each row of the new dframe and remove the old column headings
melt1$level.order <- seq(1, length(melt1$SampleID),1)
melt1 <- melt1[,c(1:6,8:9)]

# create a dframe with all the pollen counts for each sample
melt2 <- melt(dframe1[,c("SampleID","Family_Species","Date","Site","SlideNumber","PollenCount",
                         "PollenNumber_1","PollenNumber_2","PollenNumber_3","PollenNumber_4",
                         "PollenNumber_5","PollenNumber_6","PollenNumber_7","PollenNumber_8",
                         "PollenNumber_9","PollenNumber_10","PollenNumber_11","PollenNumber_12","PollenNumber_13")],
              id=c("SampleID","Family_Species","Date","Site","SlideNumber","PollenCount"))

# label each row of the new dframe and remove the old column headings
melt2$level.order <- seq(1, length(melt2$SampleID),1)
melt2 <- melt2[,c(1:6,8:9)]  #trim out the ID column - this is important for merging back together

# rename the 'value' column on the second dframe so it merges correctly
colnames(melt2) <- c("SampleID","Family_Species","Date","Site","SlideNumber","PollenCount", "value2","level.order")

# merge the dframes back together
dframe1.melt <- merge(melt1, melt2)

# remove the row labels as we are finished with them
dframe1.melt <- dframe1.melt[,c(1:6,8:9)]

# remove all the rows with NAs, keeping a single row for any moths with zero pollen (because of the zeroes in PollenNumber_1)
dframe1.melt <- dframe1.melt[complete.cases(dframe1.melt[,8]),]

# rename the remaining columns to original names
colnames(dframe1.melt) <- c("SampleID","Family_Species","Date","Site","SlideNumber","PollenCount", "PollenType","PollenNumber")

# view information about the columns
str(dframe1.melt)

# reformat the data into a matrix
matrix1 <- dcast(dframe1.melt, SampleID + Family_Species + Date +
                   Site + SlideNumber + PollenCount ~ PollenType,
                 value.var = "PollenNumber",
                 fun.aggregate = sum)



# this adds an extra blank column (I don't know why!) but let's delete it
# if you apply this to a different dataset you'll need to check the number of columns and adjust as appropriate
matrix1 <- matrix1[,c(1:6,8:78)]


### create a variable for fire/no fire.
### This tells R to search for "NF" in each row of 'Site', to label rows containing "NF" as NoFire, and all other rows as Fire
matrix1$Treatment <- ifelse(grepl("NF",matrix1$Site),"NoFire","Fire")
matrix1$Treatment <- factor(matrix1$Treatment)

### create a variable for sample
# this will be important when we come to do network analysis
# first tell R to give each date a unique value
matrix1$SamplingDay <- factor(matrix1$Date,
                              levels = levels(matrix1$Date),
                              labels = 1:21)

# then create a new column including both Site and Date information, so each sample at each site has a unique level
matrix1$Sample <- do.call(paste, c(matrix1[c("Site","SamplingDay")], sep = "_"))

# finally, reorder the columns to put pollen variables after the new variables
names(matrix1) # check what columns we have and what order
matrix1 <- matrix1[,c(1:6,78:80,7:77)] # tell R what order to put them in
names(matrix1) # check it's worked


# output the reformatted dataframe to a .txt file to use in downstream analysis
write.table(matrix1, "Data\\MatrixNoct.txt", sep="\t", row.names=FALSE)

# N.B. at this point any typing errors in species names will become obvious as they will be assigned an additional column.
# Therefore at this point, go back to the raw data, find these typing errors, correct them, and run this script again.



# now we want to produce a separate file with the data organised by each sampling session at each site
# we will do this by aggregating samples from the same site and sampling session

matrix1r <- matrix1[,c(3:4,9:80)]
summary(matrix1r)


matrix2 <- ddply(matrix1r, .(Site,Date,Sample), numcolwise(sum))

# output the reformatted dataframe to a .txt file to use in downstream analysis
write.table(matrix2, "Data\\SamplesNoct.txt", sep="\t", row.names=FALSE)

