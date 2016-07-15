###################################################################################
####   Script for basic insect abundance and assemblage composition analysis   ####
###################################################################################


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
summary(dframe3$Count)

# data appear Poisson and possibly overdispersed
# a quick check for overdispersion; in a regular Poisson distribution, mean=variance:
mean(dframe3$Count)
var(dframe3$Count)  # var > mean therefore data are overdispersed


### Data appear overdispersed, so we will try two types of model in addition to Poisson: quasi-Poisson, and negative binomial
### However not many zeroes in data, so zero-inflated models would be inappropriate

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
# you need to pass a model and you can optionally pass up to two explanatory variables as well
# you want residuals that are roughly normally distributed around zero with no obvious trends
# you want residuals that are roughly the same for each level of your treatments
chkres(model1P,dframe3$Treatment,dframe3$Season)  # these residuals are really not bad considering the potential overdispersion



# try a log-transformation
hist(log(dframe3$Count+1,10)) # this looks quite normally distributed

dframe3$tCount <- log(dframe3$Count+1,10) # log-transform the data, adding 1 (as 0 cannot be log transformed)


plot(tCount ~ Treatment, data = dframe3)
plot(tCount ~ Season, data = dframe3)
plot(tCount ~ Month, data = dframe3)



# Gaussian

model1G <- lmer(tCount ~ Treatment*Season
                + (1|Year) + (1|Site) + (1|Date), # random effects
                data = dframe3)

summary(model1G)
drop1(model1G, test="Chi")

chkres(model1G,dframe3$Treatment,dframe3$Season) # residuals are pretty good 

# but interaction term non-significant so reconstruct without this

model1Ga <- lmer(tCount ~ Treatment + Season
                 + (1|Year) + (1|Site) + (1|Date), # random effects
                 data = dframe3)

summary(model1Ga)
drop1(model1Ga, test="Chi")

chkres(model1Ga,dframe3$Treatment,dframe3$Season) # we have a winner 




### Species Richness ###

# plot richness against Treatment and Season
plot(SpeciesRichness ~ interaction(Treatment,Season), data = dframe3)


# test effects of Treatment and Season on abundance
# first plot distribution of abundance
hist(dframe3$SpeciesRichness)
hist(log(dframe3$SpeciesRichness+1, 10)) # There are zeroes in this data which can't be log transformed
summary(dframe3$SpeciesRichness)


# data are Poisson and possibly overdispersed
mean(dframe3$SpeciesRichness)
var(dframe3$SpeciesRichness)  # var > mean therefore data are overdispersed


# Poisson model using lme4

model2P <- glmer(SpeciesRichness ~ Treatment * Season # fixed effects
                 + (1|Year) + (1|Site) + (1|Date), # random effects
                 family = poisson (link = "log"),
                 data = dframe3)

# warning about model convergence so check it again, using:

chkconv(model2P) # model is acceptable

# inspect and test the model
summary(model2P)
drop1(model2P, test="Chisq")  

# check the model's residuals
chkres(model2P,dframe3$Treatment,dframe3$Season)  # these residuals aren't great for the Treatment variable

# try a QP for comparison's sake

# QuasiPoisson model using MASS

model2Q <- glmmPQL(SpeciesRichness ~ Treatment * Season,
                   random = list(~1|Year, ~1|Site, ~1|Date),
                   family = quasipoisson (link = "log"),
                   data = dframe3)


summary(model2Q)
Anova(model2Q, type="III")  # drop1 doesn't work properly with this model class so we use a Type III Anova instead

# check residuals
# this function produces a subset of the previous plots that are available for this model class
chkres.PQL(model2Q,dframe3$Treatment,dframe3$Season) # these look a little better but are positive skewed


# try a log-transformation
hist(log(dframe3$SpeciesRichness+1,10)) # this looks quite normally distributed

dframe3$tSpeciesRichness <- log(dframe3$SpeciesRichness+1,10) # log-transform the data, adding 1 (as 0 cannot be log transformed)

plot(tSpeciesRichness ~ interaction(Treatment,Season), data = dframe3)


# Gaussian

model2G <- lmer(tSpeciesRichness ~ Treatment*Season
                + (1|Year) + (1|Site) + (1|Date), # random effects
                data = dframe3)

summary(model2G)
drop1(model2G, test="Chi")

chkres(model2G,dframe3$Treatment,dframe3$Season) # residuals are fairly good 

# but interaction term marginally non-significant so reconstruct without this

model2Ga <- lmer(tSpeciesRichness ~ Treatment + Season
                 + (1|Year) + (1|Site) + (1|Date), # random effects
                 data = dframe3)#

summary(model2Ga)
drop1(model2Ga, test="Chi")

chkres(model2Ga,dframe3$Treatment,dframe3$Season) # a very slight hint of positive trend, but better than the others 





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
dframe5 <- ddply(dframe4r, .(Date,Site,Transect,Treatment,Sample,Season,Month,Year), numcolwise(mean))

dframe5SR <- ddply(dframe4r, .(Date,Site,Transect,Treatment,Sample,Season,Month,Year), numcolwise(sum))


dframe5$Month <- ordered(dframe5$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

summary(dframe5)


# Floral abundance
plot(PlantCoverage ~ Treatment, dframe5)
plot(PlantCoverage ~ Season, dframe5)
plot(PlantCoverage ~ Month, dframe5)
plot(PlantCoverage ~ Year, dframe5)
plot(PlantCoverage ~ interaction(Treatment,Season), dframe5)

hist(dframe5$PlantCoverage)
hist(log(dframe5$PlantCoverage+1,10))
summary(dframe5$PlantCoverage)

# possibly overdispersion, possibly zero-inflation (1st quartile is 0)
# because it's average data it's non-integer so Poisson won't work
# but the log looks close to Gaussian so let's try that first
mean(dframe5$PlantCoverage)
var(dframe5$PlantCoverage)  # var > mean therefore data are overdispersed

# Gaussian

dframe5$tPlantCoverage <- log(dframe5$PlantCoverage + 1,10)
summary(dframe5$tPlantCoverage)
hist(dframe5$tPlantCoverage)


model5G <- lmer(tPlantCoverage ~ Treatment*Season
                 + (1|Year) + (1|Site) + (1|Date),
                 data = dframe5)

summary(model5G)
drop1(model5G, test = "Chi")

chkres(model5G,dframe5$Treatment,dframe5$Season) # these look ok-ish 
# but have a fairly clear negative trend and problems with the fit of season


# Poisson
# can't deal with non-integers so we'll have to use total rather than mean coverage

model5P <- glmer(PlantCoverage ~ Treatment*Season
                 + (1|Year) + (1|Site) + (1|Date),
                 family = poisson (link = "log"),
                 data = dframe5SR)

summary(model5P)
drop1(model5P, test = "Chi")

chkres(model5P,dframe5SR$Treatment,dframe5SR$Season) # these are pretty bad

# nbinom

model5NB <- glmer.nb(PlantCoverage ~ Treatment*Season
                  + (1|Year) + (1|Site) + (1|Date),
                  data = dframe5SR)

chkconv(model5NB)

summary(model5NB)
drop1(model5NB, test = "Chi")

chkres(model5NB,dframe5SR$Treatment,dframe5SR$Season) # these are pretty bad as well

# we need to bite the bullet and use zero-inflated models
# zero-inflated poisson

#model5ZIP <- glmmadmb(PlantCoverage ~  Treatment*Season
#                      + (1|Year) + (1|Site) + (1|Date), #Random effects
#                      zeroInflation=TRUE,
#                      family = "poisson",
#                      link = "log",
#                      data = dframe5SR)

#summary(model5ZIP)

#Anova(model5ZIP, type="III")
#drop1(model5ZIP, test="Chisq")

#chkres.zi(model5ZIP,dframe5SR$Treatment,dframe5SR$Season) # these are pretty horrible


# zero-inflated nbinom

#model5ZINB <- glmmadmb(PlantCoverage ~  Treatment*Season
#                       + (1|Year) + (1|Site) + (1|Date), #Random effects
#                       zeroInflation=TRUE,
#                       family = "nbinom",
#                       link = "log",
#                       data = dframe5SR)

#summary(model5ZINB)

#Anova(model5ZINB, type = "III")
#drop1(model5ZINB, test="Chisq")

#chkres.zi(model5ZINB,dframe5SR$Treatment,dframe5SR$Season)



# actually everything is substantially worse than that log-Gaussian model at the start
# let's just use that one.

summary(model5G)
drop1(model5G, test = "Chi")




# Floral Species Richness
plot(SpeciesRichness ~ Treatment, dframe5SR)
plot(SpeciesRichness ~ Season, dframe5SR)
plot(SpeciesRichness ~ Month, dframe5SR)
plot(SpeciesRichness ~ interaction(Treatment,Season), dframe5SR)

hist(dframe5SR$SpeciesRichness)

# possibly overdispersion (but not sure...)
mean(dframe5SR$SpeciesRichness)
var(dframe5SR$SpeciesRichness)  # var > mean therefore data are overdispersed


# Poisson model
model6P <- glmer(SpeciesRichness ~ Treatment*Season
                 + (1|Year) + (1|Site) + (1|Date),
                 family = poisson (link="log"),
                 data = dframe5SR)

summary(model6P)
drop1(model6P, test = "Chi")

chkres(model6P,dframe5SR$Treatment,dframe5SR$Season) # these are BAAAD - they look like the data needs a log transformation


# negative binomial
model6NB <- glmer.nb(SpeciesRichness ~ Treatment*Season
                     + (1|Year) + (1|Site) + (1|Date),
                     data = dframe5SR)

summary(model6NB)
drop1(model6NB, test="Chi")

chkres(model6NB,dframe5SR$Treatment,dframe5SR$Season)


# try a log-transformed gaussian as that worked previously
dframe5SR$tSpeciesRichness <- log(dframe5SR$SpeciesRichness+1,10)
hist(dframe5SR$tSpeciesRichness)

model6G <- lmer(tSpeciesRichness ~ Treatment*Season
                + (1|Year) + (1|Site) + (1|Date),
                data = dframe5SR)

summary(model6G)
drop1(model6G, test = "Chi")

chkres(model6G,dframe5SR$Treatment,dframe5SR$Season) # these are not terrible, though the zeroes in summer continue to have a visible effect






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







