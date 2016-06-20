########################################################
####   Script for basic pollen transport analysis   ####
########################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","car","ggplot2","RVAideMemoire","arm","MASS","MuMIn","bbmle")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

# for some reason the package glmmADMB won't install via the usual methods, so:
#install.packages("R2admb")
#install.packages("glmmADMB", 
#                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                         getOption("repos")),
#                 type="source")
library(glmmADMB)



### load up Callum's custom set of functions
source("CheckResidsFunction.R")


### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/MatrixNoct.txt", header=TRUE)

summary(dframe1) # Check it's imported correctly

### tell R that SampleID and SlideNumber should be treated as factors
dframe1$SampleID <- factor(dframe1$SampleID)
dframe1$SlideNumber <- factor(dframe1$SlideNumber)

summary(dframe1)

### total up the pollen grains for each insect
dframe1$PollenLoad <- rowSums(dframe1[c(10:82)])

### total up the number of pollen species for each insect (i.e. how many columns do not contain 0)
dframe1$PollenTypes <- rowSums(dframe1[c(10:82)] != 0)

### create a binary (yes/no) variable for whether each insect is carrying any pollen
dframe1$PollenYN <- ifelse(dframe1$PollenTypes==0,0,1)



### create a subset dframe containing only the interactions
interactions <- subset(dframe1, select=-c(SampleID,Date,Site,SlideNumber,PollenCount,Treatment,SamplingDay,Sample,PollenTypes,PollenLoad,PollenYN))
summary(interactions)


### now you're ready to start looking for patterns!
  


### Let's first look at pollen load per-moth
### Plot it against treatment so you have an idea of what to expect
plot(PollenLoad ~ Treatment, data = dframe1)
hist(dframe1$PollenLoad)

### Data clearly have lots of skew, so it's worth checking for overdispersion
### Simplest test is to compare the mean and variance. 
### In a regular Poisson distribution, mean ~= variance; if variance is much larger, data are overdispersed
mean(dframe1$PollenLoad)
var(dframe1$PollenLoad)

### Data appear overdispersed, so we will try two types of model in addition to Poisson: quasi-Poisson, and negative binomial
### However no zeroes in data, so zero-inflated models would be inappropriate

# construct models using Date and Site as random effects
# random effects are factors that might affect the output variable, but not in a way that is interesting to us
# you might see variation between different sampling days, or between different fire or non-fire sites...
# ...but we are only really interested in variation due to Treatment


# Poisson model using lme4

model1P <- glmer(PollenLoad ~ Treatment # fixed effects
                 + (1|Date) + (1|Site), # random effects
                 family = poisson (link = "log"),
                 data = dframe1)

# inspect and test the model
summary(model1P)

# the best test of a GLMM is a Likelihood Ratio Test, which compares the fit of the model to if each term was removed one at a time using Chi-squared tests
# to run an LRT, use the drop1 function (i.e. drop 1 term at a time), specifying Chi-squared
drop1(model1P, test="Chisq")  



# check the model's residuals
# this custom function produces a selection of plots that you can scroll through to check that residuals look ok
# you want residuals that are roughly normally distributed around zero with no obvious trends
chkres(model1P)  # these residuals do appear to have a slight negative trend so they are not ideal


# QuasiPoisson model using MASS

model1Q <- glmmPQL(PollenLoad ~ Treatment,
                   random = list(~1|Date, ~1|Site),
                   family = quasipoisson (link = "log"),
                   data = dframe1)


summary(model1Q)
Anova(model1Q, type="III")  # drop1 doesn't work properly with this model class so we use a Type III Anova instead

# check residuals
# this function produces a subset of the previous plots that are available for this model class
chkres.PQL(model1Q) # these look slightly better but still not great



# Negative binomial model using lme4

model1NB <- glmer.nb(PollenLoad ~ Treatment # fixed effects
                     + (1|Date) + (1|Site), # random effects
                     data = dframe1)

summary(model1NB)
drop1(model1NB, test= "Chisq") # glmer.nb produces the same model class as glmer so we can treat it the same

chkres(model1NB) # these are still not ideal but they're the best yet

### choose from these candidate error families
### we can see from the residuals that model1NB is the best option so we use this result


### therefore Treatment does not significantly affect per-moth pollen load




### Let's now look at pollen types per-moth
### Plot it against treatment so you have an idea of what to expect
plot(PollenTypes ~ Treatment, data = dframe1)
hist(dframe1$PollenTypes)

### Again data clearly have lots of skew, though not so bad, so it's worth checking for overdispersion
### Simplest test is to compare the mean and variance. 
### In a regular Poisson distribution, mean ~= variance; if variance is much larger, data are overdispersed
mean(dframe1$PollenTypes)
var(dframe1$PollenTypes)

### Data appear possibly overdispersed, so we will try two types of model in addition to Poisson: quasi-Poisson, and negative binomial
### However no zeroes in data, so zero-inflated models would be inappropriate

# construct models using Date and Site as random effects


# Poisson model using lme4

model2P <- glmer(PollenTypes ~ Treatment # fixed effects
                 + (1|Date) + (1|Site), # random effects
                 family = poisson (link = "log"),
                 data = dframe1)

# inspect and test the model
summary(model2P)

# the best test of a GLMM is a Likelihood Ratio Test, which compares the fit of the model to if each term was removed one at a time using Chi-squared tests
# to run an LRT, use the drop1 function (i.e. drop 1 term at a time), specifying Chi-squared
drop1(model2P, test="Chisq")  



# check the model's residuals
chkres(model2P)  # these residuals do appear to have a very slight negative trend but probably nothing to worry about too much


# QuasiPoisson model using MASS

model2Q <- glmmPQL(PollenTypes ~ Treatment,
                   random = list(~1|Date, ~1|Site),
                   family = quasipoisson (link = "log"),
                   data = dframe1)


summary(model2Q)
Anova(model2Q, type="III")  # drop1 doesn't work properly with this model class so we use a Type III Anova instead

# check residuals
chkres.PQL(model2Q) # these look better in some ways, worse in others


# Negative binomial model using lme4

model2NB <- glmer.nb(PollenTypes ~ Treatment # fixed effects
                     + (1|Date) + (1|Site), # random effects
                     data = dframe1)

summary(model2NB)
drop1(model2NB, test= "Chisq") # glmer.nb produces the same model class as glmer so we can treat it the same

chkres(model2NB) # these are very similar to the Poisson residuals

### choose from these candidate error families
### the residuals are all similarly good. In this situation it's best to choose the simplest model - in this case, Poisson

### therefore Treatment does not significantly affect per-moth pollen types

