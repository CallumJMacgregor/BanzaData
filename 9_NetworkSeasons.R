###########################################################
####   Script for pollen transport analysis by season  ####
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



### Flower visitation ###


### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe1<-read.table("Data/NetworkMetricsFVbySeason.txt", header=TRUE)

summary(dframe1) # Check it's imported correctly

# we need a variable to allow the networks to be paired:

dframe1$ExactSeason <- do.call(paste, c(dframe1[c("Season","Year")], sep = "_"))
dframe1$ExactSeason <- as.factor(dframe1$ExactSeason)
summary(dframe1$ExactSeason)

### now you're ready to start looking for patterns!

### let's start simple - links per species

### Plot it against treatment so you have an idea of what to expect
hist(dframe1$links.per.species)                   # this looks like a Poisson distribution, but...
hist(log(dframe1$links.per.species,10))           # ... the logarithm looks like a normal distribution
plot(links.per.species ~ Treatment, data = dframe1)


# construct models using ExactSeason as random effects
# random effects are factors that might affect the output variable, but not in a way that is interesting to us
# you might see variation between different sampling days, or between different fire or non-fire sites...
# ...but we are only really interested in variation due to Treatment
# in this case, we can use the random effect to specify that our data points are paired

# although the data look Poisson, technically Poisson is only meant to handle integer values, which these are not...
# instead, let's use Gaussian distribution on the logarithm of the values

model1LPS <- lmer(log(links.per.species,10) ~ Treatment*Season # fixed effects
                  + (1|ExactSeason), # random effects
                  data = dframe1)

# inspect and test the model
summary(model1LPS)


# the best test of a GLMM is a Likelihood Ratio Test, which compares the fit of the model to if each term was removed one at a time using Chi-squared tests
# to run an LRT, use the drop1 function (i.e. drop 1 term at a time), specifying Chi-squared
drop1(model1LPS, test="Chi")  


# check the model's residuals
# this custom function produces a selection of plots that you can scroll through to check that residuals look ok
# you want residuals that are roughly normally distributed around zero (in the histogram)
# and with no obvious trends (in the scatterplots)
chkres(model1LPS,dframe1$Treatment,dframe1$Season)  # these residuals look fine


### next - linkage density

hist(dframe1$linkage.density)   # this looks pretty good for a Gaussian distribution
plot(linkage.density ~ Treatment, data = dframe1)

# construct the model

model1LD <- lmer(linkage.density ~ Treatment*Season
                 + (1|ExactSeason),
                 data=dframe1)

# inspect and test
summary(model1LD)
drop1(model1LD, test = "Chi")

chkres(model1LD,dframe1$Treatment,dframe1$Season)

# try again without non-significant interaction term


model1LDa <- lmer(linkage.density ~ Treatment+Season
                  + (1|ExactSeason),
                  data=dframe1)

# inspect and test
summary(model1LDa)
drop1(model1LDa, test = "Chi")

chkres(model1LDa,dframe1$Treatment,dframe1$Season)  # nothing concerning in these


### next - vulnerability

hist(dframe1$vulnerability.LL)  # not perfectly normal but not drastically non-normal
plot(vulnerability.LL ~ Treatment, data = dframe1)

# construct the model

model1V <- lmer(vulnerability.LL ~ Treatment*Season
                + (1|ExactSeason),
                data = dframe1)

summary(model1V)
drop1(model1V, test = "Chi")

chkres(model1V, dframe1$Treatment, dframe1$Season)


# try again without non-significant interaction term


model1Va <- lmer(vulnerability.LL ~ Treatment+Season
                  + (1|ExactSeason),
                  data=dframe1)

# inspect and test
summary(model1Va)
drop1(model1Va, test = "Chi")

chkres(model1Va,dframe1$Treatment,dframe1$Season)  # nothing concerning in these


### next - generality

hist(dframe1$generality.HL)
plot(generality.HL ~ Treatment, data = dframe1)

# construct the model

model1G <- lmer(generality.HL ~ Treatment*Season
                + (1|ExactSeason),
                data=dframe1)

summary(model1G)
drop1(model1G, test = "Chi")

chkres(model1G,dframe1$Treatment,dframe1$Season)


# try again without non-significant interaction term

model1Ga <- lmer(generality.HL ~ Treatment+Season
                 + (1|ExactSeason),
                 data=dframe1)

summary(model1Ga)
drop1(model1Ga, test = "Chi")

chkres(model1Ga,dframe1$Treatment,dframe1$Season)  # these are mostly ok



### next - interaction evenness

### Plot it against treatment so you have an idea of what to expect
hist(dframe1$interaction.evenness)                   # this looks like a normal distribution
plot(interaction.evenness ~ Treatment, data = dframe1)


# construct models using Date and Site as random effects

model1IE <- lmer(interaction.evenness ~ Treatment*Season # fixed effects
                 + (1|ExactSeason), # random effects
                 data = dframe1)

# inspect and test the model
summary(model1IE)
drop1(model1IE, test="Chi")  


# check the model's residuals
chkres(model1IE, dframe1$Treatment, dframe1$Season)  # these residuals look ok


# try again without the non-significant interaction term

model1IEa <- lmer(interaction.evenness ~ Treatment+Season # fixed effects
                  + (1|ExactSeason), # random effects
                  data = dframe1)

# inspect and test the model
summary(model1IEa)
drop1(model1IEa, test="Chi")  


# check the model's residuals
chkres(model1IEa, dframe1$Treatment, dframe1$Season)  # these residuals look ok





### next - robustness

### Plot it against treatment so you have an idea of what to expect
hist(dframe1$robustness)                   # this looks roughly like a Poisson distribution but non-integer
hist(log(dframe1$robustness,10))           # this looks roughly like a normal distribution
plot(robustness ~ Treatment, data = dframe1)


# construct models using Date and Site as random effects

model1RB <- lmer(log(robustness,10) ~ Treatment*Season # fixed effects
                 + (1|ExactSeason), # random effects
                 data = dframe1)

# inspect and test the model
summary(model1RB)
drop1(model1RB, test="Chi")  


# check the model's residuals
chkres(model1RB, dframe1$Treatment, dframe1$Season)  # these residuals aren't too bad - a hint of a positive trend


# try again without significant interaction

model1RBa <- lmer(log(robustness,10) ~ Treatment+Season
                  + (1|ExactSeason),
                  data = dframe1)

summary(model1RBa)
drop1(model1RBa, test = "Chi")

chkres(model1RBa, dframe1$Treatment, dframe1$Season)  # a hint of a positive trend and an imbalance in Treatment here



# try different error families


# inverse Gaussian (we don't often use this but conceptually it's appropriate for robustness)


model1RBig <- glmer(robustness ~ Treatment*Season # fixed effects
                    + (1|ExactSeason), # random effects
                    family = inverse.gaussian(link="log"),
                    data = dframe1)

# inspect and test the model
summary(model1RBig)
drop1(model1RBig, test="Chi")  


# check the model's residuals
chkres(model1RBig, dframe1$Treatment, dframe1$Season)  # these residuals are ok



# try quasi-poisson

model1RBqp <- glmmPQL(robustness ~ Treatment*Season,
                      random = list(~1|ExactSeason),
                      family = quasipoisson (link = "log"),
                      data = dframe1)

summary(model1RBqp)
Anova(model1RBqp, type="III")

chkres.PQL(model1RBqp, dframe1$Treatment, dframe1$Season)  # these aren't an improvement



# inverse Gaussian was the best option
drop1(model1RBig, test = "Chi")






### Pollen transport ###



### read in the data - this is the .txt file you produced in the PreparingData.R script. 
dframe2<-read.table("Data/NetworkMetricsPTbySeason.txt", header=TRUE)

summary(dframe2) # Check it's imported correctly

# we need a variable to allow the networks to be paired:

dframe2$ExactSeason <- do.call(paste, c(dframe2[c("Season","Year")], sep = "_"))
dframe2$ExactSeason <- as.factor(dframe2$ExactSeason)
summary(dframe2$ExactSeason)


### now you're ready to start looking for patterns!

### let's start simple - links per species

### Plot it against treatment so you have an idea of what to expect
hist(dframe2$links.per.species)                   # this looks like a Poisson distribution, but...
hist(log(dframe2$links.per.species,10))           # ... the logarithm looks like a normal distribution
plot(log(links.per.species,10) ~ Treatment, data = dframe2)


# construct models using Date and Site as random effects
# random effects are factors that might affect the output variable, but not in a way that is interesting to us
# you might see variation between different sampling days, or between different fire or non-fire sites...
# ...but we are only really interested in variation due to Treatment

# although the data look Poisson, technically Poisson is only meant to handle integer values, which these are not...
# instead, let's use Gaussian distribution on the logarithm of the values

model2LPS <- lmer(log(links.per.species,10) ~ Treatment*Season # fixed effects
                  + (1|ExactSeason), # random effects
                  data = dframe2)

# inspect and test the model
summary(model2LPS)


# the best test of a GLMM is a Likelihood Ratio Test, which compares the fit of the model to if each term was removed one at a time using Chi-squared tests
# to run an LRT, use the drop1 function (i.e. drop 1 term at a time), specifying Chi-squared
drop1(model2LPS, test="Chi")  


# check the model's residuals
# this custom function produces a selection of plots that you can scroll through to check that residuals look ok
# you want residuals that are roughly normally distributed around zero (in the histogram)
# and with no obvious trends (in the scatterplots)
chkres(model2LPS, dframe2$Treatment, dframe2$Season)  # these residuals look fine




### next - linkage density

hist(dframe2$linkage.density)   # this looks pretty good for a Gaussian distribution
plot(linkage.density ~ Treatment, data = dframe2)

# construct the model

model2LD <- lmer(linkage.density ~ Treatment*Season
                 + (1|ExactSeason),
                 data=dframe2)

# inspect and test
summary(model2LD)
drop1(model2LD, test = "Chi")

chkres(model2LD,dframe2$Treatment,dframe2$Season)

# try again without non-significant interaction term


model2LDa <- lmer(linkage.density ~ Treatment+Season
                  + (1|ExactSeason),
                  data=dframe2)

# inspect and test
summary(model2LDa)
drop1(model2LDa, test = "Chi")

chkres(model2LDa,dframe2$Treatment,dframe2$Season)  # nothing too worrying in these



### next - vulnerability

hist(dframe2$vulnerability.LL)  # not perfectly normal but not drastically non-normal
plot(vulnerability.LL ~ Treatment, data = dframe2)

# construct the model

model2V <- lmer(vulnerability.LL ~ Treatment*Season
                + (1|ExactSeason),
                data = dframe2)

summary(model2V)
drop1(model2V, test = "Chi")

chkres(model2V, dframe2$Treatment, dframe2$Season)


# try again without non-significant interaction term


model2Va <- lmer(vulnerability.LL ~ Treatment+Season
                 + (1|ExactSeason),
                 data=dframe2)

# inspect and test
summary(model2Va)
drop1(model2Va, test = "Chi")

chkres(model2Va,dframe2$Treatment,dframe2$Season)  # nothing concerning in these


### next - generality

hist(dframe2$generality.HL)
plot(generality.HL ~ Treatment, data = dframe2)

# construct the model

model2G <- lmer(generality.HL ~ Treatment*Season
                + (1|ExactSeason),
                data=dframe2)

summary(model2G)
drop1(model2G, test = "Chi")

chkres(model2G,dframe2$Treatment,dframe2$Season)


# try again without non-significant interaction term

model2Ga <- lmer(generality.HL ~ Treatment+Season
                 + (1|ExactSeason),
                 data=dframe2)

summary(model2Ga)
drop1(model2Ga, test = "Chi")

chkres(model2Ga,dframe2$Treatment,dframe2$Season)  # these are mostly ok



### next - interaction evenness

### Plot it against treatment so you have an idea of what to expect
hist(dframe2$interaction.evenness)                   # this looks like a normal distribution
plot(interaction.evenness ~ Treatment, data = dframe2)


# construct models using Date and Site as random effects

model2IE <- lmer(interaction.evenness ~ Treatment*Season # fixed effects
                 + (1|ExactSeason), # random effects
                 data = dframe2)

# inspect and test the model
summary(model2IE)
drop1(model2IE, test="Chi")  


# check the model's residuals
chkres(model2IE, dframe2$Treatment, dframe2$Season)  # these residuals look ok


# try again without the non-significant interaction term

model2IEa <- lmer(interaction.evenness ~ Treatment+Season # fixed effects
                  + (1|ExactSeason), # random effects
                  data = dframe2)

# inspect and test the model
summary(model2IEa)
drop1(model2IEa, test="Chi")  


# check the model's residuals
chkres(model2IEa, dframe2$Treatment, dframe2$Season)  # these residuals look ok





### next - robustness

### Plot it against treatment so you have an idea of what to expect
hist(dframe2$robustness)                   # this looks roughly like a Poisson distribution but non-integer
hist(log(dframe2$robustness,10))           # this looks roughly like a normal distribution
plot(robustness ~ Treatment, data = dframe2)


# construct models using Date and Site as random effects

model2RB <- lmer(log(robustness,10) ~ Treatment*Season # fixed effects
                 + (1|ExactSeason), # random effects
                 data = dframe2)

# inspect and test the model
summary(model2RB)
drop1(model2RB, test="Chi")  


# check the model's residuals
chkres(model2RB, dframe2$Treatment, dframe2$Season)  # these residuals aren't too bad - a hint of a positive trend



# try different error families


# inverse Gaussian (we don't often use this but conceptually it's appropriate for robustness)


model2RBig <- glmer(robustness ~ Treatment*Season # fixed effects
                    + (1|ExactSeason), # random effects
                    family = inverse.gaussian(link="log"),
                    data = dframe2)

# inspect and test the model
summary(model2RBig)
drop1(model2RBig, test="Chi")  


# check the model's residuals
chkres(model2RBig, dframe2$Treatment, dframe2$Season)  # these residuals are ok - they show less imbalance in residual variance between treatments



# try quasi-poisson

model2RBqp <- glmmPQL(robustness ~ Treatment*Season,
                      random = list(~1|ExactSeason),
                      family = quasipoisson (link = "log"),
                      data = dframe2)

summary(model2RBqp)
Anova(model2RBqp, type="III")

chkres.PQL(model2RBqp, dframe2$Treatment, dframe2$Season)  # these aren't an improvement



# inverse Gaussian was the best option
drop1(model2RBig, test = "Chi")




