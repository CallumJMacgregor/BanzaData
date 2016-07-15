###############################################################################
####   Script for analysis of functionally extrapolated species richness   ####
###############################################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","car","ggplot2","RVAideMemoire","arm","MASS","plyr","reshape2","vegan")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

# for some reason the package glmmADMB won't install via the usual methods, so:
#install.packages("R2admb")
#install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")),type="source")
library(glmmADMB)



### load up Callum's custom set of functions
f <- c("CheckResidsFunction.R","CheckConvergenceFunction.R","SpecAccumFunctions.R")
lapply(f, source)


### Moths

### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/MatrixNoct.txt", header=TRUE)

### tell R that SampleID and SamplingDay should be treated as factors
dframe1$SamplingDay <- factor(dframe1$SamplingDay)
dframe1$Month <- ordered(dframe1$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
dframe1$Year <- factor(dframe1$Year)

summary(dframe1) # Check it's imported correctly


# We don't need the pollen data or the individual SampleID for this analysis, so get rid of it:

dframe1r <- dframe1[,c(2:4,7:12)]
summary(dframe1r)

# each row currently contains 1 insect, so add a column for "Count"; 
# this is 1 in every instance except the one sample with zero insects
dframe1r$Count <- ifelse(dframe1r$Family_Species=="none",0,1)

# summarise the dataframe to how many individuals of each species in each sample
dframe2 <- ddply(dframe1r, .(Family_Species,Site,Treatment,Date,Sample,Season,Month,Year), numcolwise(sum))



### Functional extrapolation of species richness


# we have two options - either estimate the site-based richness based on occurrance in multiple samples,
# or estimate the sample-based richness based on abundance in the sample

# for the moth data, the latter is probably better due to high seasonal turnover
# however, we could do the former for site_season-based richness as well - ie accumulation for each site across years but within seasons
# doing both has the benefit of letting us see if we get roughly equivalent results

dframe2$SiteSeason <- do.call(paste, c(dframe2[c("Site","Season")], sep = "_"))
dframe2$SiteSeason <- as.factor(dframe2$SiteSeason)

dframe2$TreatmentSeason <- do.call(paste, c(dframe2[c("Treatment","Season","Year")], sep = "_"))
dframe2$TreatmentSeason <- as.factor(dframe2$TreatmentSeason)

# dframe2 contains information of how many of each species are in each sample
# we'll need this information in a matrix format though


matrix1 <- dcast(dframe2, Date + Site + Treatment
                 + Sample + Season + Month + Year + SiteSeason + TreatmentSeason
                 ~ Family_Species,
                 value.var = "Count",
                 fun.aggregate = sum)

# we need to remove the column for "none", which is the 242nd
matrix1 <- matrix1[,c(1:241,243:length(matrix1))]


# first, let's inspect the species accumulation curves for the data
rownames(matrix1) <- matrix1[,4]
matrix1a <- matrix1[,-c(1:9)]

# ...accumulation for each site within each season
matrices1 <- split(matrix1, list(matrix1$SiteSeason))  # this creates a list of smaller dframes, one for each level of sample

# first let's use another of Callum's custom functions to plot all the species accumulation curves
lapply(matrices1, sacplot)

# let's also try for each site irrespective of season
matrices2 <- split(matrix1, list(matrix1$Site))  # this creates a list of smaller dframes, one for each level of sample
lapply(matrices2, sacplot)



# we can see that very few of these show any signs of nearing an asymptote, even with larger numbers of samples

# therefore we need to extrapolate species richness
# first try sample-based - doing it for every sample based on abundance within the sample
matrices3 <- split(matrix1, list(matrix1$Sample))
SampleSR <- lapply(matrices3, samplebased)

### we now have a list of dataframes, each one containing the sample-level SR of one sample

# for sites that had too little data, the dataframes contain the value "NA" for one of the estimators
# combine all dataframes together
SampleSR.merge <- do.call("cbind", SampleSR)    # merge the data with one sample per column
colnames(SampleSR.merge) <- names(SampleSR)     # assign the sample names to each column
SampleSR.merge <- data.frame(t(SampleSR.merge))
SampleSR.merge


### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the sample-level pollen dataframe:
dframeP<-read.table("Data/SamplesNoct.txt", header=TRUE)
summary(dframeP)

# we don't need any of the pollen data, so:
dframeP <- dframeP[,1:8]

# but we do need a Treatment variable, so:
dframeP$Treatment <- ifelse(grepl("NF",dframeP$Site),"NoFire","Fire")
dframeP$Treatment <- factor(dframeP$Treatment)

# for merge to work, we need to set the row names
rownames(dframeP) <- dframeP$Sample

# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
SampleSR.full <- merge(dframeP, SampleSR.merge, by=0)





# now try site-based - doing it for every site based on repeated sampling at the site
SiteSR <- lapply(matrices1, sitebased)

### we now have a list of dataframes, each one containing the sample-level SR of one sample

# for sites that had too little data, the dataframes contain the value "NA" for one of the estimators
# combine all dataframes together
SiteSR.merge <- do.call("cbind", SiteSR)    # merge the data with one sample per column
colnames(SiteSR.merge) <- names(SiteSR)     # assign the sample names to each column
SiteSR.merge <- data.frame(t(SiteSR.merge))
SiteSR.merge


### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the sample-level pollen dataframe:
# but it needs collapsing to the same number of sites
dframePS <- matrix1[,c(1:10)]
dframePS <- ddply(dframePS, .(Treatment,Season,SiteSeason), numcolwise(sum))
dframePS <- dframePS[,c(1:3)]

# for merge to work, we need to set the row names
rownames(dframePS) <- dframePS$SiteSeason

# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
SiteSR.full <- merge(dframePS, SiteSR.merge, by=0)




# finally try by Treatment + Season (18 paired samples)
matrices4 <- split(matrix1, list(matrix1$TreatmentSeason))
TreatmentSR <- lapply(matrices4, sitebased)

### we now have a list of dataframes, each one containing the sample-level SR of one sample

# for sites that had too little data, the dataframes contain the value "NA" for one of the estimators
# combine all dataframes together
TreatmentSR.merge <- do.call("cbind", TreatmentSR)    # merge the data with one sample per column
colnames(TreatmentSR.merge) <- names(TreatmentSR)     # assign the sample names to each column
TreatmentSR.merge <- data.frame(t(TreatmentSR.merge))

### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is held within the sample-level pollen dataframe:
# but it needs collapsing to the same number of sites
dframePST <- matrix1[,c(1:10)]
dframePST <- ddply(dframePST, .(Treatment,Season,TreatmentSeason), numcolwise(sum))
dframePST <- dframePST[,c(1:3)]

# for merge to work, we need to set the row names
rownames(dframePST) <- dframePST$TreatmentSeason

# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
TreatmentSR.full <- merge(dframePST, TreatmentSR.merge, by=0)




### let's analyse each of these in turn, using the Chao1 estimated SR

# first by sample
summary(SampleSR.full)

hist(SampleSR.full$S.chao1)
hist(log(SampleSR.full$S.chao1+1, 10))

plot(S.chao1 ~ Treatment, SampleSR.full)
plot(S.chao1 ~ Season, SampleSR.full)
plot(S.chao1 ~ interaction(Treatment,Season), SampleSR.full)


model1G <- lmer(log(S.chao1+1,10) ~ Treatment + Season
                + (1|Year) + (1|Site) + (1|Date),
                data = SampleSR.full)

summary(model1G)
drop1(model1G, test = "Chi")

chkres(model1G,SampleSR.full$Treatment,SampleSR.full$Season)  # these are actually not at all bad



# by site
summary(SiteSR.full)

plot(chao ~ Treatment, SiteSR.full)
plot(chao ~ Season, SiteSR.full)
plot(chao ~ interaction(Treatment,Season), SiteSR.full)


hist(SiteSR.full$chao)
hist(log(SiteSR.full$chao+1,10))

model2QP <- glm(chao ~ Treatment + Season,
                family = quasipoisson (link = "log"),
                data=SiteSR.full)

summary(model2QP)
drop1(model2QP, test = "Chi")

chkres(model2QP,SiteSR.full$Treatment,SiteSR.full$Season)




# by treatment
summary(TreatmentSR.full)

# in order for this to be treated as a paired test we need to create a factor with a unique level for each exact season

# first break the TreatmentSeason variable up at the _ breaks
TreatmentSR.full <- separate(data = TreatmentSR.full, col = TreatmentSeason, into = c("Treatment1", "ExactSeason","Year"), sep = "_")

# then stitch the Season and Year components back together
TreatmentSR.full$ExactSeason <- do.call(paste, c(TreatmentSR.full[c("ExactSeason","Year")], sep = "_"))
TreatmentSR.full$ExactSeason <- as.factor(TreatmentSR.full$ExactSeason)

# get rid of the columns we don't need
TreatmentSR.full <- TreatmentSR.full[,c(1:3,5,7:length(TreatmentSR.full))]

plot(chao ~ Treatment, TreatmentSR.full)
plot(chao ~ Season, TreatmentSR.full)
plot(chao ~ interaction(Treatment,Season), TreatmentSR.full)


hist(TreatmentSR.full$chao) # this definitely looks like Poisson - but values are non-integer
# it's also not horribly different to a Gaussian

# Gaussian
# random effect of the ExactSeason variable lets model treat this as a paired test

model3G <- lmer(chao ~ Treatment*Season
                + (1|ExactSeason),
                data=TreatmentSR.full)

summary(model3G)
drop1(model3G, test = "Chi")

chkres(model3G,TreatmentSR.full$Treatment,TreatmentSR.full$Season) # some imbalance in Season

# interaction non-sig so let's try without


model3Ga <- lmer(chao ~ Treatment + Season
                 +(1|ExactSeason),
                 data=TreatmentSR.full)

summary(model3Ga)
drop1(model3Ga, test = "Chi")

chkres(model3Ga,TreatmentSR.full$Treatment,TreatmentSR.full$Season)
# slightly better though still some imbalance in Season, also residuals are not normally distributed


# try a quasipoisson

model3QP <- glmmPQL(chao ~ Treatment + Season,
                    random = ~1|ExactSeason,
                    family = quasipoisson (link = "log"),
                    data=TreatmentSR.full)

summary(model3QP)
Anova(model3QP, type = "III")

chkres.PQL(model3QP,TreatmentSR.full$Treatment,TreatmentSR.full$Season)

# still would be nice to try it as a Poisson, so let's try rounding SR estimates to the nearest integer

TreatmentSR.full$chao.int <- round(TreatmentSR.full$chao, digits = 0)

plot(chao.int ~ Treatment, TreatmentSR.full)
plot(chao.int ~ Season, TreatmentSR.full)
plot(chao.int ~ interaction(Treatment,Season), TreatmentSR.full)


hist(TreatmentSR.full$chao.int) # this definitely looks like Poisson - and values are now integer
mean(TreatmentSR.full$chao.int)
var(TreatmentSR.full$chao.int)

model3P <- glmer(chao.int ~ Treatment * Season
                 +(1|ExactSeason),
                 family = poisson (link = "log"),
                 data = TreatmentSR.full)

summary(model3P)
drop1(model3P, test = "Chi")

chkres(model3P,TreatmentSR.full$Treatment,TreatmentSR.full$Season) # these aren't terrible but not perfect either


# try a neg binom as this looks a tad overdispersed

model3NB <- glmer.nb(chao.int ~ Treatment + Season
                     +(1|ExactSeason),
                     data = TreatmentSR.full)

chkconv(model3NB)

summary(model3NB)
drop1(model3NB, test = "Chi")

chkres(model3NB,TreatmentSR.full$Treatment,TreatmentSR.full$Season)

# Poisson is the most appropriate model in theory and the fit is acceptable, so we proceed with it.
summary(model3P)
drop1(model3P, test = "Chi")



### Plants

### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe3<-read.table("Data/PlantTransects.txt", header=TRUE)

dframe3$Year <- factor(dframe3$Year)
summary(dframe3)


# We don't need the individual TransectID for this analysis, so get rid of it:

dframe3r <- dframe3[,c(2:14)]
summary(dframe3r)

# each row currently contains 1 plant species, so add a column for "Count"; 
# this is 1 in every instance except the transect with zero insects
dframe3r$Transects <- ifelse(dframe3r$PlantSpecies=="none",0,1)

# summarise the dataframe to how much total coverage of each species in each sample (non-integers could be a problem so we won't use mean)
dframe4 <- ddply(dframe3r, .(Date,Site,Treatment,Sample,Season,Month,Year,PlantSpecies), numcolwise(sum))

# set up the SiteSeason and TreatmentSeason variables for condensing the data later

dframe4$SiteSeason <- do.call(paste, c(dframe4[c("Site","Season")], sep = "_"))
dframe4$SiteSeason <- as.factor(dframe4$SiteSeason)


dframe4$TreatmentSeason <- do.call(paste, c(dframe4[c("Treatment","Season","Year")], sep = "_"))
dframe4$TreatmentSeason <- as.factor(dframe4$TreatmentSeason)

# put the data out into a matrix format

matrix2 <- dcast(dframe4, Date + Site + Treatment
                 + Sample + Season + Month + Year + SiteSeason + TreatmentSeason
                 ~ PlantSpecies,
                 value.var = "PlantCoverage",
                 fun.aggregate = sum)

# we need to remove the column for "none", which is the 52nd and the mystery blank, which is the 10th
matrix2 <- matrix2[,c(1:9,11:51,53:length(matrix2))]

# first, let's inspect the species accumulation curves for the data
rownames(matrix2) <- matrix2[,4]
matrix2a <- matrix2[,-c(1:9)]

# ...accumulation for each site within each season
matricesP1 <- split(matrix2, list(matrix2$SiteSeason))  # this creates a list of smaller dframes, one for each level of sample


# first let's use another of Callum's custom functions to plot all the species accumulation curves
lapply(matricesP1, sacplot)

# let's also try for each site irrespective of season
matricesP2 <- split(matrix2, list(matrix2$Site))  # this creates a list of smaller dframes, one for each level of sample
lapply(matricesP2, sacplot)


# some of these might be nearing an asymptote but none are there, even with larger numbers of samples

# therefore we need to extrapolate species richness

# sample-based won't work here as it's dependent on the number of singletons
# as this is percent cover, 1% cover would be a singleton, which isn't really appropriate


# now try site-based - doing it for every site based on repeated sampling at the site
SitePSR <- lapply(matricesP1, sitebased)

### we now have a list of dataframes, each one containing the sample-level SR of one sample

# for sites that had too little data, the dataframes contain the value "NA" for one of the estimators
# combine all dataframes together
SitePSR.merge <- do.call("cbind", SitePSR)    # merge the data with one sample per column
colnames(SitePSR.merge) <- names(SitePSR)     # assign the sample names to each column
SitePSR.merge <- data.frame(t(SitePSR.merge))

### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is previously prepared in dframePS:


# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
SitePSR.full <- merge(dframePS, SitePSR.merge, by=0)



# finally try by Treatment + Season (18 paired samples)
matricesP3 <- split(matrix2, list(matrix2$TreatmentSeason))
TreatmentPSR <- lapply(matricesP3, sitebased)

# combine all dataframes together
TreatmentPSR.merge <- do.call("cbind", TreatmentPSR)    # merge the data with one sample per column
colnames(TreatmentPSR.merge) <- names(TreatmentPSR)     # assign the sample names to each column
TreatmentPSR.merge <- data.frame(t(TreatmentPSR.merge))


### before this can be analysed, we need to reintroduce data about each site (e.g. Treatment)
### all this information is previously prepared in dframePST:

# now we merge the two dataframes together by row name, using 'by=0' (column 0 is the row names)
TreatmentPSR.full <- merge(dframePST, TreatmentPSR.merge, by=0)




### let's analyse each of these in turn, using the Chao1 estimated SR

# first, sites

summary(SitePSR.full)

plot(chao ~ Treatment, SitePSR.full)
plot(chao ~ Season, SitePSR.full)
plot(chao ~ interaction(Treatment,Season), SitePSR.full)


hist(SitePSR.full$chao)
hist(log(SitePSR.full$chao,10))


model4QP <- glm(chao ~ Treatment + Season,
                family = quasipoisson (link = "log"),
                data=SitePSR.full)

summary(model4QP)
drop1(model4QP, test = "Chi")

chkres(model4QP,SitePSR.full$Treatment,SitePSR.full$Season)  # really decent fit here





# second, paired treatment-seasons

summary(TreatmentPSR.full)

# in order for this to be treated as a paired test we need to create a factor with a unique level for each exact season

# first break the TreatmentSeason variable up at the _ breaks
TreatmentPSR.full <- separate(data = TreatmentPSR.full, col = TreatmentSeason, into = c("Treatment1", "ExactSeason","Year"), sep = "_")

# then stitch the Season and Year components back together
TreatmentPSR.full$ExactSeason <- do.call(paste, c(TreatmentPSR.full[c("ExactSeason","Year")], sep = "_"))
TreatmentPSR.full$ExactSeason <- as.factor(TreatmentPSR.full$ExactSeason)

# get rid of the columns we don't need
TreatmentPSR.full <- TreatmentPSR.full[,c(1:3,5,7:length(TreatmentPSR.full))]



plot(chao ~ Treatment, TreatmentPSR.full)
plot(chao ~ Season, TreatmentPSR.full)
plot(chao ~ interaction(Treatment,Season), TreatmentPSR.full)


hist(TreatmentPSR.full$chao)
hist(log(TreatmentPSR.full$chao,10))

model5QP <- glmmPQL(chao ~ Treatment + Season,
                    random = ~1|ExactSeason,
                    family = quasipoisson (link = "log"),
                    data=TreatmentPSR.full)

summary(model5QP)
Anova(model5QP, type = "III")

chkres.PQL(model5QP,TreatmentPSR.full$Treatment,TreatmentPSR.full$Season)  # not tooo bad!



model5G <- lmer(log(chao,10) ~ Treatment * Season
                + (1|ExactSeason),
                data=TreatmentPSR.full)

summary(model5G)
drop1(model5G, test = "Chi")

chkres(model5G,TreatmentSR.full$Treatment,TreatmentSR.full$Season) # not bad


# try it as a Poisson, so let's try rounding SR estimates to the nearest integer

TreatmentPSR.full$chao.int <- round(TreatmentPSR.full$chao, digits = 0)

plot(chao.int ~ Treatment, TreatmentPSR.full)
plot(chao.int ~ Season, TreatmentPSR.full)
plot(chao.int ~ interaction(Treatment,Season), TreatmentPSR.full)


hist(TreatmentPSR.full$chao.int) # this definitely looks like Poisson - and values are now integer
mean(TreatmentPSR.full$chao.int)
var(TreatmentPSR.full$chao.int)

model5P <- glmer(chao.int ~ Treatment * Season
                 +(1|ExactSeason),
                 family = poisson (link = "log"),
                 data = TreatmentPSR.full)

summary(model5P)
drop1(model5P, test = "Chi")

chkres(model5P,TreatmentPSR.full$Treatment,TreatmentPSR.full$Season) # these aren't terrible but not perfect either









### Community dissimilarity using Adonis (vegan)



### Moths

# we need the data in the original matrix format
summary(matrix1)

rownames(matrix1) <- matrix1$Sample

# but let's trim off some of the excess columns - we only need Treatment, Season and Site
matrix1d <- matrix1[,c(2:3,5,10:length(matrix1))]

# we need to get rid of any rows with zeroes in this time
matrix1d$Total <- rowSums(matrix1d[,c(4:length(matrix1d))])

summary(matrix1d$Total)

matrix1ds <- subset(matrix1d,matrix1d$Total>0)

# we also need a version with all factor columns removed
community1 <- matrix1ds[,c(4:330)]

# first we want to visualise things so as to know what to expect
mod1 <- metaMDS(community1)

MothsTreat <- plot(mod1)
ordihull(mod1, group=matrix1ds$Treatment, show="Fire")
ordihull(mod1, group=matrix1ds$Treatment, show="NoFire")

MothsSeason <- plot(mod1)
ordihull(mod1, group=matrix1ds$Season, show="Spring")
ordihull(mod1, group=matrix1ds$Season, show="Summer")
ordihull(mod1, group=matrix1ds$Season, show="Autumn")
ordihull(mod1, group=matrix1ds$Season, show="Winter")

<<<<<<< HEAD

=======
>>>>>>> a9f0f32237963c4baf4315babbf86dce32c2d398
# based on these plots there may be an effect of Treatment and there almost certainly will be an effect of Season

# now we use adonis to test the dissimilarities between communities in different samples
# we will ask it to test for effects of Treatment and Season (and an interaction), constraining permutations to within Sites
# we'll set the method as Bray-Curtis

AdMoths <- adonis(community1 ~ Treatment*Season,
                  data = matrix1ds,
                  strata = matrix1ds$Site,
                  method = "bray",
                  perm=1e5)

AdMoths

# interaction non-significant but repeating without doesn't really change things




### Plants


# we need the data in the original matrix format
summary(matrix2)


# but let's trim off some of the excess columns - we only need Treatment, Season and Site
matrix2d <- matrix2[,c(2:3,5,10:length(matrix2))]





# we need to get rid of any rows with zeroes in this time
matrix2d$Total <- rowSums(matrix2d[,c(4:length(matrix2d))])

summary(matrix2d$Total)

matrix2ds <- subset(matrix2d,matrix2d$Total>0)





# we also need a version with all factor columns removed
community2 <- matrix2ds[,c(4:74)]
community2[is.na(community2)] <- 0


# first we want to visualise things so as to know what to expect
mod2 <- metaMDS(community2)

FlowersTreat <- plot(mod2)
ordihull(mod2, group=matrix2ds$Treatment, show="Fire")
ordihull(mod2, group=matrix2ds$Treatment, show="NoFire")

MothsSeason <- plot(mod2)
ordihull(mod2, group=matrix2ds$Season, show="Spring")
ordihull(mod2, group=matrix2ds$Season, show="Summer")
ordihull(mod2, group=matrix2ds$Season, show="Autumn")
ordihull(mod2, group=matrix2ds$Season, show="Winter")

# these plots aren't that informative!

# now we use adonis to test the dissimilarities between communities in different samples
# we will ask it to test for effects of Treatment and Season (and an interaction), constraining permutations to within Sites
# we'll set the method as Bray-Curtis

AdFlowers <- adonis(community2 ~ Treatment*Season,
                  data = matrix2ds,
                  strata = matrix2ds$Site,
                  method = "bray",
                  perm=1e5)

AdFlowers






############################ development ##################################



