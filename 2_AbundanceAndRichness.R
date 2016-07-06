###################################################################################
####   Script for basic insect abundance and assemblage composition analysis   ####
###################################################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","car","ggplot2","RVAideMemoire","arm","MASS","plyr")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

# for some reason the package glmmADMB won't install via the usual methods, so:
#install.packages("R2admb")
#install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")),type="source")
library(glmmADMB)



### load up Callum's custom set of functions
source("CheckResidsFunction.R")


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

# and additionally summarise it into how many insects (of any species) in each sample
dframe2r <- dframe2[,-1]
dframe2r$SpeciesRichness <- ifelse(dframe2r$Count==0,0,1)
dframe3 <- ddply(dframe2r, .(Site,Treatment,Date,Sample,Season,Month,Year), numcolwise(sum))


# now we're ready to start analysing 
summary(dframe3)

# plot abundance against Treatment and Season
plot(Count ~ Treatment, data = dframe3)
plot(Count ~ Season, data = dframe3)
plot(Count ~ Month, data = dframe3)

# test effects of Treatment and Season on abundance
# first plot distribution of abundance
hist(dframe3$Count)

# data are Poisson and possibly overdispersed
# a quick check for overdispersion; in a regular Poisson distribution, mean=variance:
mean(dframe3$Count)
var(dframe3$Count)  # var > mean therefore data are overdispersed


### Data appear overdispersed, so we will try two types of model in addition to Poisson: quasi-Poisson, and negative binomial
### However no zeroes in data, so zero-inflated models would be inappropriate

# construct models using Date and Site as random effects
# random effects are factors that might affect the output variable, but not in a way that is interesting to us
# you might see variation between different sampling days, or between different fire or non-fire sites...
# ...but we are only really interested in variation due to Treatment


# Poisson model using lme4

model1P <- glmer(Count ~ Treatment * Season # fixed effects
                 + (1|Year) + (1|Site) + (1|Date), # random effects
                 family = poisson (link = "log"),
                 data = dframe3)

# inspect and test the model
summary(model1P)

# the best test of a GLMM is a Likelihood Ratio Test, which compares the fit of the model to if each term was removed one at a time using Chi-squared tests
# to run an LRT, use the drop1 function (i.e. drop 1 term at a time), specifying Chi-squared
drop1(model1P, test="Chisq")  



# check the model's residuals
# this custom function produces a selection of plots that you can scroll through to check that residuals look ok
# you want residuals that are roughly normally distributed around zero with no obvious trends
chkres(model1P)  # these residuals are really not bad considering the potential overdispersion


# QuasiPoisson model using MASS

model1Q <- glmmPQL(Count ~ Treatment * Season,
                   random = list(~1|Year, ~1|Site, ~1|Date),
                   family = quasipoisson (link = "log"),
                   data = dframe3)


summary(model1Q)
Anova(model1Q, type="III")  # drop1 doesn't work properly with this model class so we use a Type III Anova instead

# check residuals
# this function produces a subset of the previous plots that are available for this model class
chkres.PQL(model1Q) # these look markedly worse



# Negative binomial model using lme4

model1NB <- glmer.nb(Count ~ Treatment * Season # fixed effects
                     + (1|Year) + (1|Site) + (1|Date), # random effects
                     data = dframe3)

summary(model1NB)
drop1(model1NB, test= "Chisq") # glmer.nb produces the same model class as glmer so we can treat it the same

chkres(model1NB) # these are worse again

### choose from these candidate error families
### none are perfect so we use the simplest model with best fit - Poisson

summary(model1P)
drop1(model1P, test="Chisq")


### Species Richness ###

# plot richness against Treatment and Season
plot(SpeciesRichness ~ Treatment, data = dframe3)
plot(SpeciesRichness ~ Season, data = dframe3)


# test effects of Treatment and Season on abundance
# first plot distribution of abundance
hist(dframe3$SpeciesRichness)

# data are Poisson and possibly overdispersed
# a quick check for overdispersion; in a regular Poisson distribution, mean=variance:
mean(dframe3$SpeciesRichness)
var(dframe3$SpeciesRichness)  # var > mean therefore data are overdispersed


# Poisson model using lme4

model2P <- glmer(SpeciesRichness ~ Treatment * Season # fixed effects
                 + (1|Year) + (1|Site) + (1|Date), # random effects
                 family = poisson (link = "log"),
                 data = dframe3)

# inspect and test the model
summary(model2P)

# the best test of a GLMM is a Likelihood Ratio Test, which compares the fit of the model to if each term was removed one at a time using Chi-squared tests
# to run an LRT, use the drop1 function (i.e. drop 1 term at a time), specifying Chi-squared
drop1(model2P, test="Chisq")  



# check the model's residuals
# this custom function produces a selection of plots that you can scroll through to check that residuals look ok
# you want residuals that are roughly normally distributed around zero with no obvious trends
chkres(model2P)  # these residuals do appear to have a slight negative trend so they are not ideal



### Plants

### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe4<-read.table("Data/PlantTransects.txt", header=TRUE)

dframe4$Year <- factor(dframe4$Year)
summary(dframe4)


# We don't need the individual TransectID for this analysis, so get rid of it:

dframe4r <- dframe4[,c(2:14)]
summary(dframe4r)

# each row currently contains 1 plant species, so add a column for "Count"; 
# this is 1 in every instance except the transect with zero insects
dframe4r$SpeciesRichness <- ifelse(dframe4r$PlantSpecies=="none",0,1)

# summarise the dataframe to how many individuals of each species in each sample
dframe5 <- ddply(dframe4r, .(Date,Site,Transect,Treatment,Sample,Season,Month,Year), numcolwise(sum))

dframe5$Month <- ordered(dframe5$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

summary(dframe5)


# Floral abundance
plot(PlantCoverage ~ Treatment, dframe5)
plot(PlantCoverage ~ Season, dframe5)
plot(PlantCoverage ~ Month, dframe5)
plot(PlantCoverage ~ interaction(Treatment,Season), dframe5)

hist(dframe5$PlantCoverage)

# possibly overdispersion (but not sure...)
mean(dframe5$PlantCoverage)
var(dframe5$PlantCoverage)  # var > mean therefore data are overdispersed


# Poisson model
model5P <- glmer(PlantCoverage ~ Treatment*Season
                 + (1|Year) + (1|Site) + (1|Date),
                 family = poisson (link="log"),
                 data = dframe5)

summary(model5P)
drop1(model5P, test = "Chi")

chkres(model5P)


# Floral Species Richness
plot(SpeciesRichness ~ Treatment, dframe5)
plot(SpeciesRichness ~ Season, dframe5)
plot(SpeciesRichness ~ Month, dframe5)
plot(SpeciesRichness ~ interaction(Treatment,Season), dframe5)

hist(dframe5$SpeciesRichness)

# possibly overdispersion (but not sure...)
mean(dframe5$SpeciesRichness)
var(dframe5$SpeciesRichness)  # var > mean therefore data are overdispersed


# Poisson model
model6P <- glmer(SpeciesRichness ~ Treatment*Season
                 + (1|Year) + (1|Site) + (1|Date),
                 family = poisson (link="log"),
                 data = dframe5)

summary(model6P)
drop1(model6P, test = "Chi")

chkres(model6P)









##################development###################

dframe3Sum <- dframe3[dframe3$Season=="Summer",]
dframe3Spr <- dframe3[dframe3$Season=="Spring",]
dframe3Win <- dframe3[dframe3$Season=="Winter",]

plot(Count ~ Treatment, dframe3Sum)
plot(Count ~ Treatment, dframe3Spr)
plot(Count ~ Treatment, dframe3Win)


modelWinP  <- glmer(Count ~ Treatment # fixed effects
                    + (1|Site), # random effects
                    family = poisson (link = "log"),
                    data = dframe3Win)
summary(modelWinP)
drop1(modelWinP, test="Chi")


plot(SpeciesRichness ~ Treatment, dframe3Sum)
plot(SpeciesRichness ~ Treatment, dframe3Spr)
plot(SpeciesRichness ~ Treatment, dframe3Win)
