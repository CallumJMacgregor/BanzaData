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

### Insects and pollen

### read in the data - this is the raw data as input by Paula. 
# I have manually added a column called 'SampleID' with a unique number for each row of the dataframe;
# I have renamed "Family/Species" to "Family_Species", 
# and added _1 to the first PollenType and PollenNumber column for consistency
# Finally, I've put a 0 in "PollenNumber_1" for insects with no pollen - this is important
dframe1<-read.csv("Data/BanzaDataNoct.csv", header=TRUE)

summary(dframe1) # Check it's imported correctly

### Remove redundant PollenTypes 14 & 15
dframe1 <- dframe1[,1:33]
summary(dframe1)

### Make binary "Lepidoptera" variable a factor
dframe1$Lepidoptera <- factor(dframe1$Lepidoptera)

### Melt the dataframe
names(dframe1) # check what columns you have and what they're called

# create a dframe with all the family/species names for each sample
melt1 <- melt(dframe1[,c("SampleID","Family_Species","Lepidoptera","Date","Site","SlideNumber","PollenCount",
                         "PollenType_1","PollenType_2","PollenType_3","PollenType_4",
                         "PollenType_5","PollenType_6","PollenType_7","PollenType_8",
                         "PollenType_9","PollenType_10","PollenType_11","PollenType_12","PollenType_13")],
                      id=c("SampleID","Family_Species","Lepidoptera","Date","Site","SlideNumber","PollenCount"))

#label each row of the new dframe and remove the old column headings
melt1$level.order <- seq(1, length(melt1$SampleID),1)
melt1 <- melt1[,-8]

# create a dframe with all the pollen counts for each sample
melt2 <- melt(dframe1[,c("SampleID","Family_Species","Lepidoptera","Date","Site","SlideNumber","PollenCount",
                         "PollenNumber_1","PollenNumber_2","PollenNumber_3","PollenNumber_4",
                         "PollenNumber_5","PollenNumber_6","PollenNumber_7","PollenNumber_8",
                         "PollenNumber_9","PollenNumber_10","PollenNumber_11","PollenNumber_12","PollenNumber_13")],
              id=c("SampleID","Family_Species","Lepidoptera","Date","Site","SlideNumber","PollenCount"))

# label each row of the new dframe and remove the old column headings
melt2$level.order <- seq(1, length(melt2$SampleID),1)
melt2 <- melt2[,-8]  #trim out the ID column - this is important for merging back together

# rename the 'value' column on the second dframe so it merges correctly
colnames(melt2) <- c("SampleID","Family_Species","Lepidoptera","Date","Site","SlideNumber","PollenCount", "value2","level.order")

# merge the dframes back together
dframe1.melt <- merge(melt1, melt2)

# remove the row labels as we are finished with them
dframe1.melt <- dframe1.melt[,-8]

# remove all the rows with NAs, keeping a single row for any moths with zero pollen (because of the zeroes in PollenNumber_1)
dframe1.melt <- dframe1.melt[complete.cases(dframe1.melt[,9]),]

# rename the remaining columns to original names
colnames(dframe1.melt) <- c("SampleID","Family_Species","Lepidoptera","Date","Site","SlideNumber","PollenCount", "PollenType","PollenNumber")

# view information about the columns
str(dframe1.melt)

# reformat the data into a matrix
matrix1 <- dcast(dframe1.melt, SampleID + Family_Species + Lepidoptera + Date +
                   Site + SlideNumber + PollenCount ~ PollenType,
                 value.var = "PollenNumber",
                 fun.aggregate = sum)



# this adds an extra blank column (I don't know why!) and a column for zero; let's delete them now
matrix1 <- matrix1[,c(1:7,10:length(matrix1))]


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


# create a variable for season and one for month
# four seasons of three months each:
# "winter" - Dec-Feb, "spring" - Mar-May, "Summer" - Jun-Aug, "Autumn" - Sep-Nov

# this code searches in the date column for the month and assigns each month to the appropriate categories

matrix1$Season <- ifelse(grepl("/11/",matrix1$Date),"Autumn",
                         ifelse(grepl("/12/",matrix1$Date),"Autumn",
                                ifelse(grepl("/01/",matrix1$Date),"Winter",
                                       ifelse(grepl("/02/",matrix1$Date),"Winter",
                  ifelse(grepl("/03/",matrix1$Date),"Winter",
                         ifelse(grepl("/04/",matrix1$Date),"Spring",
                                ifelse(grepl("/05/",matrix1$Date),"Spring",
                                       ifelse(grepl("/06/",matrix1$Date),"Spring",
                  ifelse(grepl("/07/",matrix1$Date),"Summer",
                         ifelse(grepl("/08/",matrix1$Date),"Summer",
                                ifelse(grepl("/09/",matrix1$Date),"Summer",
                                       ifelse(grepl("/10/",matrix1$Date),"Autumn",
                                              "Fail"))))))))))))
matrix1$Season <- factor(matrix1$Season)

matrix1$Month <- ifelse(grepl("/11/",matrix1$Date),"Nov",
                         ifelse(grepl("/12/",matrix1$Date),"Dec",
                                ifelse(grepl("/01/",matrix1$Date),"Jan",
                                       ifelse(grepl("/02/",matrix1$Date),"Feb",
                                              ifelse(grepl("/03/",matrix1$Date),"Mar",
                                                     ifelse(grepl("/04/",matrix1$Date),"Apr",
                                                            ifelse(grepl("/05/",matrix1$Date),"May",
                                                                   ifelse(grepl("/06/",matrix1$Date),"Jun",
                                                                          ifelse(grepl("/07/",matrix1$Date),"Jul",
                                                                                 ifelse(grepl("/08/",matrix1$Date),"Aug",
                                                                                        ifelse(grepl("/09/",matrix1$Date),"Sep",
                                                                                               ifelse(grepl("/10/",matrix1$Date),"Oct",
                                                                                                      "Fail"))))))))))))
matrix1$Month<-ordered(matrix1$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))


matrix1$Year <- ifelse(grepl("/04/2013",matrix1$Date),"1",
                       ifelse(grepl("/05/2013",matrix1$Date),"1",
                              ifelse(grepl("/06/2013",matrix1$Date),"1",
                                     ifelse(grepl("/07/2013",matrix1$Date),"1",
                                            ifelse(grepl("/08/2013",matrix1$Date),"1",
                ifelse(grepl("/09/2014",matrix1$Date),"3",
                       ifelse(grepl("/10/2014",matrix1$Date),"3",
                              ifelse(grepl("/11/2014",matrix1$Date),"3",
                                     ifelse(grepl("/12/2014",matrix1$Date),"3",
                                            ifelse(grepl("/01/2015",matrix1$Date),"3",
                                                   ifelse(grepl("/02/2015",matrix1$Date),"3",
                                                          ifelse(grepl("/03/2015",matrix1$Date),"3",
                                                                 ifelse(grepl("/04/2015",matrix1$Date),"3",
                                                                        ifelse(grepl("/05/2015",matrix1$Date),"3",
                                                    "2"))))))))))))))
matrix1$Year <- factor(matrix1$Year)


# check it's worked - there should be no values for Fail
summary(matrix1$Season)
summary(matrix1$Month)
summary(matrix1$Year)


# finally, reorder the columns to put pollen variables after the new variables
names(matrix1) # check what columns we have and what order
matrix1 <- matrix1[,c(1:7,79:84,8:78)] # tell R what order to put them in
names(matrix1) # check it's worked


# we want to check and remove "environmental contamination" - i.e. Pinus, Briza and Cupressus
# First, to check the quantity of it:

# create a reduced dataframe
matrix1c <- matrix1[,c(4:5,10:length(matrix1))]
summary(matrix1c)

# create a dummy variable that should be identical for every row of data
matrix1c$dummy <- ifelse(grepl("Fail",matrix1c$Season),"Fail","Pass")
matrix1c$dummy <- factor(matrix1c$dummy)
summary(matrix1c$dummy)


# create a dataframe with the total quantity of each pollen type
matrix1ce <- ddply(matrix1c, .(dummy), numcolwise(sum))
summary(matrix1ce)

# calculate the percentage of total pollen grains that are (a) Pinus and (b) Cupressus
matrix1ce$Total <- rowSums(matrix1ce[,2:length(matrix1ce)])

matrix1ce$percPinus <- matrix1ce$`Pinus sp.`*100/matrix1ce$Total
matrix1ce$percCupressus <- matrix1ce$`Cupressus sp.`*100/matrix1ce$Total
matrix1ce$percBriza <- matrix1ce$`Briza maxima`*100/matrix1ce$Total

# view percentages
str(matrix1ce$percPinus)
str(matrix1ce$percCupressus)
str(matrix1ce$percBriza)


# now remove these columns from the main working dataframe
matrix1$`Pinus sp.` <- NULL
matrix1$`Cupressus sp.` <- NULL
matrix1$`Briza maxima` <- NULL


# finally, we want to restrict the dataset to Lepidoptera only
matrix1L <- matrix1[matrix1$Lepidoptera==1,]
matrix1L <- matrix1L[,-3]   

# output the reformatted dataframe to a .txt file to use in downstream analysis
write.table(matrix1L, "Data\\MatrixNoct.txt", sep="\t", row.names=FALSE)

# N.B. at this point any typing errors in species names will become obvious as they will be assigned an additional column.
# Therefore at this point, go back to the raw data, find these typing errors, correct them, and run this script again.



# now we want to produce a separate file with the data organised by each sampling session at each site
# we will do this by aggregating samples from the same site and sampling session

matrix1r <- matrix1L[,c(3:4,7:length(matrix1L))]
summary(matrix1r)


matrix2 <- ddply(matrix1r, .(Site,Date,Treatment,SamplingDay,Sample,Season,Month,Year), numcolwise(sum))

# output the reformatted dataframe to a .txt file to use in downstream analysis
write.table(matrix2, "Data\\SamplesNoct.txt", sep="\t", row.names=FALSE)






### plants

### read in the data - this is the raw data as input by Paula. 
# I have manually added a column called 'TransectID' with a unique number for each row of the dataframe;
# I have renamed some of the columns, 
# and added _1 to the first Plant and PlantCoverage column for consistency
# Finally, I've put a 0 in "PlantCoverage_1" for transects with no flowers - this is important
dframe2<-read.csv("Data/PlantTransects.csv", header=TRUE)

summary(dframe2) # Check it's imported correctly

### Remove redundant Plant columns 9 & 10
dframe2 <- dframe2[,1:22]
summary(dframe2)


### Melt the dataframe
names(dframe2) # check what columns you have and what they're called

# create a dframe with all the family/species names for each sample
melt3 <- melt(dframe2[,c("TransectID","Date","Site","Transect","PercentCoverage","AverageHeight_cm",
                         "Plant_1","Plant_2","Plant_3","Plant_4",
                         "Plant_5","Plant_6","Plant_7","Plant_8")],
              id=c("TransectID","Date","Site","Transect","PercentCoverage","AverageHeight_cm"))

#label each row of the new dframe and remove the old column headings
melt3$level.order <- seq(1, length(melt3$TransectID),1)
melt3 <- melt3[,-7]

# create a dframe with all the pollen counts for each sample
melt4 <- melt(dframe2[,c("TransectID","Date","Site","Transect","PercentCoverage","AverageHeight_cm",
                         "PlantCoverage_1","PlantCoverage_2","PlantCoverage_3","PlantCoverage_4",
                         "PlantCoverage_5","PlantCoverage_6","PlantCoverage_7","PlantCoverage_8")],
              id=c("TransectID","Date","Site","Transect","PercentCoverage","AverageHeight_cm"))

# label each row of the new dframe and remove the old column headings
melt4$level.order <- seq(1, length(melt4$TransectID),1)
melt4 <- melt4[,-7]  #trim out the ID column - this is important for merging back together

# rename the 'value' column on the second dframe so it merges correctly
colnames(melt4) <- c("TransectID","Date","Site","Transect","PercentCoverage","AverageHeight_cm","value2","level.order")

# merge the dframes back together
dframe2.melt <- merge(melt3, melt4)

# remove the row labels as we are finished with them
dframe2.melt <- dframe2.melt[,-7]

# remove all the rows with NAs, keeping a single row for any moths with zero pollen (because of the zeroes in PollenNumber_1)
dframe2.melt <- dframe2.melt[complete.cases(dframe2.melt[,8]),]

# rename the remaining columns to original names
colnames(dframe2.melt) <- c("TransectID","Date","Site","Transect","PercentCoverage","AverageHeight_cm", "PlantSpecies","PlantCoverage")

# view information about the columns
str(dframe2.melt)

# reformat the data into a matrix
matrix3 <- dframe2.melt



### create a variable for fire/no fire.
### This tells R to search for "NF" in each row of 'Site', to label rows containing "NF" as NoFire, and all other rows as Fire
matrix3$Treatment <- ifelse(grepl("NF",matrix3$Site),"NoFire","Fire")
matrix3$Treatment <- factor(matrix3$Treatment)

### create a variable for sample
# this will be important when we come to do network analysis
# first tell R to give each date a unique value
matrix3$SamplingDay <- factor(matrix3$Date,
                              levels = levels(matrix3$Date),
                              labels = 1:21)

# then create a new column including both Site and Date information, so each sample at each site has a unique level
matrix3$Sample <- do.call(paste, c(matrix3[c("Site","SamplingDay")], sep = "_"))


# create a variable for season and one for month
# four seasons of three months each:
# "winter" - Dec-Feb, "spring" - Mar-May, "Summer" - Jun-Aug, "Autumn" - Sep-Nov

# this code searches in the date column for the month and assigns each month to the appropriate categories

matrix3$Season <- ifelse(grepl("/11/",matrix3$Date),"Autumn",
                         ifelse(grepl("/12/",matrix3$Date),"Autumn",
                                ifelse(grepl("/01/",matrix3$Date),"Winter",
                                       ifelse(grepl("/02/",matrix3$Date),"Winter",
                                              ifelse(grepl("/03/",matrix3$Date),"Winter",
                                                     ifelse(grepl("/04/",matrix3$Date),"Spring",
                                                            ifelse(grepl("/05/",matrix3$Date),"Spring",
                                                                   ifelse(grepl("/06/",matrix3$Date),"Spring",
                                                                          ifelse(grepl("/07/",matrix3$Date),"Summer",
                                                                                 ifelse(grepl("/08/",matrix3$Date),"Summer",
                                                                                        ifelse(grepl("/09/",matrix3$Date),"Summer",
                                                                                               ifelse(grepl("/10/",matrix3$Date),"Autumn",
                                                                                                      "Fail"))))))))))))
matrix3$Season <- factor(matrix3$Season)

matrix3$Month <- ifelse(grepl("/11/",matrix3$Date),"Nov",
                        ifelse(grepl("/12/",matrix3$Date),"Dec",
                               ifelse(grepl("/01/",matrix3$Date),"Jan",
                                      ifelse(grepl("/02/",matrix3$Date),"Feb",
                                             ifelse(grepl("/03/",matrix3$Date),"Mar",
                                                    ifelse(grepl("/04/",matrix3$Date),"Apr",
                                                           ifelse(grepl("/05/",matrix3$Date),"May",
                                                                  ifelse(grepl("/06/",matrix3$Date),"Jun",
                                                                         ifelse(grepl("/07/",matrix3$Date),"Jul",
                                                                                ifelse(grepl("/08/",matrix3$Date),"Aug",
                                                                                       ifelse(grepl("/09/",matrix3$Date),"Sep",
                                                                                              ifelse(grepl("/10/",matrix3$Date),"Oct",
                                                                                                     "Fail"))))))))))))
matrix3$Month<-ordered(matrix3$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))


matrix3$Year <- ifelse(grepl("/04/2013",matrix3$Date),"1",
                       ifelse(grepl("/05/2013",matrix3$Date),"1",
                              ifelse(grepl("/06/2013",matrix3$Date),"1",
                                     ifelse(grepl("/07/2013",matrix3$Date),"1",
                                            ifelse(grepl("/08/2013",matrix3$Date),"1",
                                                   ifelse(grepl("/09/2014",matrix3$Date),"3",
                                                          ifelse(grepl("/10/2014",matrix3$Date),"3",
                                                                 ifelse(grepl("/11/2014",matrix3$Date),"3",
                                                                        ifelse(grepl("/12/2014",matrix3$Date),"3",
                                                                               ifelse(grepl("/01/2015",matrix3$Date),"3",
                                                                                      ifelse(grepl("/02/2015",matrix3$Date),"3",
                                                                                             ifelse(grepl("/03/2015",matrix3$Date),"3",
                                                                                                    ifelse(grepl("/04/2015",matrix3$Date),"3",
                                                                                                           ifelse(grepl("/05/2015",matrix3$Date),"3",
                                                                                                                  "2"))))))))))))))
matrix3$Year <- factor(matrix3$Year)


# check it's worked - there should be no values for Fail
summary(matrix3$Season)
summary(matrix3$Month)
summary(matrix3$Year)

matrix3$PlantSpecies <- ifelse(matrix3$PlantCoverage==0,"none",matrix3$PlantSpecies)


# finally, reorder the columns to put pollen variables after the new variables
names(matrix3) # check what columns we have and what order
matrix3 <- matrix3[,c(1:6,9:14,7:8)] # tell R what order to put them in
names(matrix3) # check it's worked



# output the reformatted dataframe to a .txt file to use in downstream analysis
write.table(matrix3, "Data\\PlantTransects.txt", sep="\t", row.names=FALSE)

# N.B. at this point any typing errors in species names will become obvious as they will be assigned an additional column.
# Therefore at this point, go back to the raw data, find these typing errors, correct them, and run this script again.







#################################### development #########################


matrix1cet <- matrix1ce[,c(2:72)]

matrix1cet <- t(matrix1cet)  

d <- matrix1cet
names <- rownames(d)
rownames(d) <- NULL
data <- cbind(names,d)
data <- data.frame(data)

colnames(data) <- c("Species","Count")

data$CountN <- as.numeric(as.character(data$Count))

plot(log(CountN,10) ~ Species, data = data)
hist(data$Count)

summary(matrix1cet)
