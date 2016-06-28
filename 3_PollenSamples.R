###########################################################
####   Script for pollen transport analysis by sample  ####
###########################################################


### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","car","ggplot2","RVAideMemoire","arm","MASS","MuMIn","glmmADMB")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

# for some reason the package glmmADMB won't install via the usual methods, so if it's not already installed:
#install.packages("R2admb")
#install.packages("glmmADMB", 
#                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                         getOption("repos")),
#                 type="source")



### load up Callum's custom set of functions
source("CheckResidsFunction.R")


### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/SamplesNoct.txt", header=TRUE)

summary(dframe1) # Check it's imported correctly


### total up the pollen grains for each sample
dframe1$PollenCount <- rowSums(dframe1[c(4:74)])

### total up the number of pollen species for each sample (i.e. how many columns do not contain 0)
dframe1$PollenTypes <- rowSums(dframe1[c(4:74)] != 0)


### create a variable for fire/no fire.
### This tells R to search for "NF" in each row of 'Site', to label rows containing "NF" as NoFire, and all other rows as Fire
dframe1$Treatment <- ifelse(grepl("NF",dframe1$Site),"NoFire","Fire")
dframe1$Treatment <- factor(dframe1$Treatment)



### now you're ready to start looking for patterns!



### Let's first look at pollen load per-moth
### Plot it against treatment so you have an idea of what to expect
plot(PollenCount ~ Treatment, data = dframe1)
hist(dframe1$PollenCount)

### Data clearly have lots of skew, so it's worth checking for overdispersion
### Simplest test is to compare the mean and variance. 
### In a regular Poisson distribution, mean ~= variance; if variance is much larger, data are overdispersed
mean(dframe1$PollenCount)
var(dframe1$PollenCount)

### Data appear overdispersed, so we will try two types of model in addition to Poisson: quasi-Poisson, and negative binomial
### However no zeroes in data, so zero-inflated models would be inappropriate

# construct models using Date and Site as random effects
# random effects are factors that might affect the output variable, but not in a way that is interesting to us
# you might see variation between different sampling days, or between different fire or non-fire sites...
# ...but we are only really interested in variation due to Treatment


# Poisson model using lme4

model1P <- glmer(PollenCount ~ Treatment # fixed effects
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
chkres(model1P)  # surprisingly despite the overdispersion these residuals look fine


# QuasiPoisson model using MASS

model1Q <- glmmPQL(PollenCount ~ Treatment,
                   random = list(~1|Date, ~1|Site),
                   family = quasipoisson (link = "log"),
                   data = dframe1)


summary(model1Q)
Anova(model1Q, type="III")  # drop1 doesn't work properly with this model class so we use a Type III Anova instead

# check residuals
# this function produces a subset of the previous plots that are available for this model class
chkres.PQL(model1Q) # these look no better, no worse



# Negative binomial model using lme4

model1NB <- glmer.nb(PollenCount ~ Treatment # fixed effects
                     + (1|Date) + (1|Site), # random effects
                     data = dframe1)

summary(model1NB)
drop1(model1NB, test= "Chisq") # glmer.nb produces the same model class as glmer so we can treat it the same

chkres(model1NB) # these are probably the worst yet

### choose from these candidate error families
### we can see from the residuals that model1P is the surprisingly the best option so we use this result


### therefore Treatment does significantly affects per-sample pollen count




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
chkres(model2P)  # These look pretty good


# QuasiPoisson model using MASS

model2Q <- glmmPQL(PollenTypes ~ Treatment,
                   random = list(~1|Date, ~1|Site),
                   family = quasipoisson (link = "log"),
                   data = dframe1)


summary(model2Q)
Anova(model2Q, type="III")  # drop1 doesn't work properly with this model class so we use a Type III Anova instead

# check residuals
chkres.PQL(model2Q) # these also look fine


# Negative binomial model using lme4

model2NB <- glmer.nb(PollenTypes ~ Treatment # fixed effects
                     + (1|Date) + (1|Site), # random effects
                     data = dframe1)

summary(model2NB)
drop1(model2NB, test= "Chisq") # glmer.nb produces the same model class as glmer so we can treat it the same

chkres(model2NB) # these are the worst of the three with indications of a trend

### choose from these candidate error families
### the residuals for Poisson and Quasi-Poisson are similarly good.
### In this situation it's best to choose the simplest model - in this case, Poisson

### therefore Treatment does significantly affect per-sample pollen types
