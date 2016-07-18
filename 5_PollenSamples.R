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
k <- c("CheckResidsFunction.R","CheckConvergenceFunction.R")
lapply(k,source)

### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/SamplesNoct.txt", header=TRUE)

summary(dframe1) # Check it's imported correctly

dframe1$Month<-ordered(dframe1$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

### total up the pollen grains for each sample
dframe1$PollenCount <- rowSums(dframe1[c(9:76)])

### total up the number of pollen species for each sample (i.e. how many columns do not contain 0)
dframe1$PollenTypes <- rowSums(dframe1[c(9:76)] != 0)



### now you're ready to start looking for patterns!



### Let's first look at pollen load per-moth
### Plot it against treatment so you have an idea of what to expect
plot(PollenCount ~ Treatment, data = dframe1)
plot(PollenCount ~ Season, data = dframe1)
hist(dframe1$PollenCount)
hist(log(dframe1$PollenCount+1,10))

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

model1P <- glmer(PollenCount ~ Treatment*Season # fixed effects
                 + (1|Year) + (1|Site) + (1|Date), # random effects
                 family = poisson (link = "log"),
                 data = dframe1)

# inspect and test the model
summary(model1P)

# the best test of a GLMM is a Likelihood Ratio Test, which compares the fit of the model to if each term was removed one at a time using Chi-squared tests
# to run an LRT, use the drop1 function (i.e. drop 1 term at a time), specifying Chi-squared
drop1(model1P, test="Chi")  



# check the model's residuals
# this custom function produces a selection of plots that you can scroll through to check that residuals look ok
# you want residuals that are roughly normally distributed around zero with no obvious trends
chkres(model1P, dframe1$Treatment, dframe1$Season)  # surprisingly despite the overdispersion these residuals look ok-ish


# QuasiPoisson model using MASS

model1Q <- glmmPQL(PollenCount ~ Treatment * Season,
                   random = list(~1|Year, ~1|Site, ~1|Date),
                   family = quasipoisson (link = "log"),
                   data = dframe1)


summary(model1Q)
Anova(model1Q, type="III")  # drop1 doesn't work properly with this model class so we use a Type III Anova instead

# check residuals
# this function produces a subset of the previous plots that are available for this model class
chkres.PQL(model1Q, dframe1$Treatment, dframe1$Season) # these look worse



# Negative binomial model using lme4

model1NB <- glmer.nb(PollenCount ~ Treatment * Season # fixed effects
                     + (1|Year) + (1|Site) + (1|Date), # random effects
                     data = dframe1)

chkconv(model1NB)

summary(model1NB)
drop1(model1NB, test= "Chi") # glmer.nb produces the same model class as glmer so we can treat it the same

chkres(model1NB, dframe1$Treatment, dframe1$Season) # these are almost ok but there does seem to be a positive trend at low values


# Gaussian with log transformation

model1G <- lmer(log(PollenCount+1,10) ~ Treatment * Season
                + (1|Year) + (1|Site) + (1|Date),
                data = dframe1)

summary(model1G)
drop1(model1G, test = "Chi")

chkres(model1G, dframe1$Treatment, dframe1$Season) # these are the best so far, but interaction is non-significant so try without


model1Ga <- lmer(log(PollenCount+1,10) ~ Treatment + Season
                 + (1|Year) + (1|Site) + (1|Date),
                 data = dframe1)

summary(model1Ga)
drop1(model1Ga, test = "Chi")

chkres(model1Ga, dframe1$Treatment, dframe1$Season) # these are the best so far




### choose from these candidate error families
### we can see from the residuals that model1Ga is the best option




### Let's now look at pollen types per-moth
### Plot it against treatment so you have an idea of what to expect
plot(PollenTypes ~ Treatment, data = dframe1)
plot(PollenTypes ~ Season, data = dframe1)
plot(PollenTypes ~ interaction(Treatment,Season), data = dframe1)
hist(dframe1$PollenTypes)
hist(log(dframe1$PollenTypes+1,10))

### Again data clearly have lots of skew, though it does look quite like a standard poisson, so it's worth checking for overdispersion
### Simplest test is to compare the mean and variance. 
### In a regular Poisson distribution, mean ~= variance; if variance is much larger, data are overdispersed
mean(dframe1$PollenTypes)
var(dframe1$PollenTypes)

# construct models using Year and Site and Date as random effects


# Poisson model using lme4

model2P <- glmer(PollenTypes ~ Treatment*Season # fixed effects
                 + (1|Year) + (1|Site) + (1|Date), # random effects
                 family = poisson (link = "log"),
                 data = dframe1)

# inspect and test the model
summary(model2P)

# the best test of a GLMM is a Likelihood Ratio Test, which compares the fit of the model to if each term was removed one at a time using Chi-squared tests
# to run an LRT, use the drop1 function (i.e. drop 1 term at a time), specifying Chi-squared
drop1(model2P, test="Chi")  



# check the model's residuals
chkres(model2P, dframe1$Treatment, dframe1$Season)  # These look pretty good but there's a hint of skew


# QuasiPoisson model using MASS

model2Q <- glmmPQL(PollenTypes ~ Treatment * Season,
                   random = list(~1|Year, ~1|Site, ~1|Date),
                   family = quasipoisson (link = "log"),
                   data = dframe1)


summary(model2Q)
Anova(model2Q, type="III")  # drop1 doesn't work properly with this model class so we use a Type III Anova instead

# check residuals
chkres.PQL(model2Q, dframe1$Treatment, dframe1$Season) # these are pretty good


# Negative binomial model using lme4

model2NB <- glmer.nb(PollenTypes ~ Treatment * Season # fixed effects
                     + (1|Year) + (1|Site) + (1|Date), # random effects
                     data = dframe1)

summary(model2NB)
drop1(model2NB, test= "Chisq") # glmer.nb produces the same model class as glmer so we can treat it the same

chkres(model2NB, dframe1$Treatment, dframe1$Season) # these are the worst of the three with indications of a trend


# Gaussian with log transformation

model2G <- lmer(log(PollenTypes+1,10) ~ Treatment * Season
                + (1|Year) + (1|Site) + (1|Date),
                data = dframe1)

summary(model2G)
drop1(model2G, test = "Chi")

chkres(model2G, dframe1$Treatment, dframe1$Season)


### choose from these candidate error families
### the residuals for Gaussian are slightly better than those for quasipoisson therefore use model2G


















############################ development #########################


dframe1c <- ddply(dframe1, .(Season,Month), numcolwise(sum))

plot(PollenCount ~ Month, data = dframe1c)
