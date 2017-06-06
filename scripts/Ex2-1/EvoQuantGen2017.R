
# This R script was written by Patrick Carter in June 2017 and is a modification of an R Scrpt written in April 2017 for Biology 521 and which itself is a modification of an R script
# written in August 2016 for the NIMBioS Evolutionary Genetics course.
# The purpose of this script is to highlight code needed to run quantitative genetic analyses in R using MCMCglmm
# The data to be analyzed are 3 phenotypic traits from a toy data set: Ptype1, Ptype2 and Ptype3.
# Data were generated for only one generation but pedigree information is known for parents and grandparents
# Individuals were measured in one of two batches

#IMPORTANT NOTE: You will need to change the path for opening data files to whatever is appropriate to your computer


#######################################################
#Visualize and examine the data 

#Clear memory if needed
rm(list=ls()) 
setwd("~/Projects/UW_EQG_2017/data/Ex2-1")

#load graphics library ggplot2
library(ggplot2)

#open the data file
Toy4 <- read.table("Toy4SimC.dat",header = T)

# look at distributions of the 3 Ptypes
Toy41.hist = hist(Toy4$Ptype1, breaks = 20)
Toy42.hist = hist(Toy4$Ptype2, breaks = 20)
Toy43.hist = hist(Toy4$Ptype3, breaks = 20)

#plot the phenotypes with each other
qplot(x = Ptype1, y = Ptype2, data = Toy4, geom = "point", color = factor(Batch) )
qplot(x = Ptype1, y = Ptype3, data = Toy4, geom = "point", color = factor(Batch) )
qplot(x = Ptype2, y = Ptype3, data = Toy4, geom = "point", color = factor(Batch) )

#Look at covariance structure and correlations
#Choose variables for which we want covariances
TempVar <- subset(Toy4, select = c(Ptype1,Ptype2,Ptype3))
summary(TempVar)
#Get covariance matrix  
CV = cov(TempVar)
#Get correlation matrix
CR = cor(TempVar)
#Look at eigen structure just for fun
EigCV <- eigen(CV)
EigCV


#######################################################
# create pedigree file 

# Important information about the structure and construction of pedigree files to be used by MCMCglmm
# pedigree files contain identification information for measured individuals and their mothers (Dam) and fathers (Sire) and other ancestors
# this file may contain many generations of individuals for whom phenotypes were never measured
# note that offspring always must have larger id numbers than their parents, parents must have larger id numbers than grandparents, etc

# make sure missing values for Sire and Dam identification numbers are coded as NA
# make sure that dam and sire id numbers are coded as factors
# the identification variable for the individuals with phenotypic information in the data set MUST be called animal and it should remain a numeric variable
# the pedigree file should be sorted by the variable animal

# individuals may appear in the pedigree as animal and they may also appear as a Sire or Dam

# frequently you must create the pedigree file from the data file, as we have to do here
# the data file contains one record for each measured individual
# each record contains data on all 3 phenotypes as well as identification information on parents and grandparents

#the MasterBayes and plyr libraries are needed to make the pedigree
library(MasterBayes)
library(plyr)

#open the data file
Toy4 <- read.table ("Toy4SimD.dat",header = T)

#Creat Offspring pedigree
Toy4Pedid<-subset(Toy4,select=c(Id,Sireid,Damid))
Toy4Pedid<-rename(Toy4Pedid,c(Id="animal"))

#Creat Pedigree of dams
#save one record per dam
Toy4PedDam<- aggregate(Toy4, list(Toy4$Damid), FUN=head, 1)
#choose only id information for pedigree
Toy4PedDam<-subset(Toy4PedDam,select=c(Damid,MGsireid,MGdamid))
#rename variables
Toy4PedDam<-rename(Toy4PedDam,c(Damid="animal",MGsireid="Sireid",MGdamid="Damid"))

#Create Pedigree of sires
#save one record per sire
Toy4PedSire<- aggregate(Toy4, list(Toy4$Sireid), FUN=head, 1)
#choose only id information for pedigree
Toy4PedSire<-subset(Toy4PedSire,select=c(Sireid,PGsireid,PGdamid))
#rename variables
Toy4PedSire<-rename(Toy4PedSire,c(Sireid="animal",PGsireid="Sireid",PGdamid="Damid"))

#Combine all 3 files using rbind to make one full pedigree file
Toy4ped<-rbind(Toy4Pedid,Toy4PedDam,Toy4PedSire)
#sort by id number
Toy4ped<-Toy4ped[order(Toy4ped$animal),]

# this commmand from the MasterBayes library completes a pedigree with missing information for some sires and dams 
# by adding the generation in which all Dams and Sires were unknown; you will need to do this when making your pedigree file:
Toy4ped<-insertPed(Toy4ped, founders=NULL) #this fills in sires and dams for ALL animals starting from 1 (parent/ grandparent generation)

#replace "." with NA for missing values for Sire and Dam id numbers
Toy4ped$Sireid[Toy4ped$Sireid=="."]<-NA
Toy4ped$Damid[Toy4ped$Damid=="."]<-NA

#identify Dam and Sire as factors
Toy4ped$Damid <- as.factor(Toy4ped$Damid)
Toy4ped$Sireid <- as.factor(Toy4ped$Sireid)

#save pedigree file as an R data file
save(Toy4ped,file="Toy4ped.rat")



###############################################
# modify data file

# make sure missing values for Sire and Dam id numbers are coded as NA
# make sure that dam and sire id numbers are coded as factors
# the identifiation variable for individuals with phenotypic data MUST called animal and it should remain a numeric variable
# drop cases of individuals without phenotypic data 
# file should be sorted by the variable animal

Toy4 <- read.table ("Toy4SimD.dat",header = T)  #dont forget to change the path#

#replace "." with NA for missing values for Sire and Dam id numbers
Toy4$Sireid[Toy4$Sireid=="."]<-NA
Toy4$Damid[Toy4$Damid=="."]<-NA

#identify Dam and Sire as factors
Toy4$Damid <- as.factor(Toy4$Damid)
Toy4$Sireid <- as.factor(Toy4$Sireid)

#rename the variable id as animal
Toy4<-rename(Toy4,c(Id="animal"))

#extract only the variables we need
Toy4<-subset(Toy4, select = c(animal,Sireid,Damid,Batch, Ptype1, Ptype2, Ptype3, Ftype1, Ftype2, Ftype3))

#summarize data file (Note: before this point you should already have thoroughly graphed and examined your data)
summary(Toy4)
#save data file as an R data file
save(Toy4,file="Toy4.rat")

#######################################################
# Model 101: Trait = Ptype1, Batch fixed, Additive effects only random effect

#needed library for running genearlized linear mixed models
library(MCMCglmm)

#Clear memory if needed
rm(list=ls()) 

#load the data file and the pedigree file
load("Toy4.rat")
load("Toy4ped.rat")

# set up parameters for the priors
# Fixed effects are automatically set by MCMCglmm to follow a normal distribution and do not need to be specified
# G is for specified random effects (additive, maternal, etc). 
# Here we set for weak priors that will be used in an inverse gamma disribution automatically set by MCMCglmm
# R is residual effects for each specified random effect and follows same rules as G
prior101 <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, n=0.002)))

#Set up inverse of Numerator Relationship Matrix using MCMCglmm inverseA command
invA <- inverseA(Toy4ped) 

# model statement
model101 <- MCMCglmm(Ftype1 ~ 1 + Batch,                #intercept (the 1) and Batch are the fixed effect
                     random = ~animal,                  #additive effects (animal) the only random effect
                     family = "gaussian",               #phenotype has gaussian distribution
                     prior = prior101,                  #call the priors parameters defined above
                     data = Toy4,                       #call the data file
                     nitt = 100000,                     #number of MCMC iterations: NORMALLY NEED ~10,000,000 for a good run and convergence. here is only 100,000 (1/100 optimal, for time's sake)
                     burnin = 10000,                    #number of iterations for burnin
                     thin = 4000,                       #sampling interval
                     ginverse=list(animal=invA)$invA)   #call the inverse of the NRM

#save model output as an R object so we can access it later; this is very important when the nitt is very high
save(model101, file = "model101.obj") #dont forget to change the path#

#load model file
load(file = "model101.obj") #dont forget to change the path#

#plot trace and density of fixed effects; should be no trend in trace
plot(model101$Sol)
#plot trace and density of random (additive and residual (=environmental)) variances; should be no trend in trace
plot(model101$VCV)

#examine autocorrelation of fixed effects
autocorr.diag(model101$Sol)
#examine autocorrelation of random (additive and residual) variances
autocorr.diag(model101$VCV)

#check effective population size for fixed effects; should be gt 1000
effectiveSize(model101$Sol)
#check effective population size for random effects (additve and residual variances); should be gt 1000
effectiveSize(model101$VCV)

#test of convergence, p should be greater than 0.05 for good convergence
heidel.diag(model101$VCV)

#estimates of additive and residual variances
posterior.mode(model101$VCV)

#summary of model; make sure to check DIC score (smaller is better)
summary(model101)

#estimate posterior distribution of the heritability (animal variance divided by animal + residual variances)
herit <- model101$VCV[, "animal"]/(model101$VCV[, "animal"] + model101$VCV[, "units"])

#effective sample size for heritability should be gt 1000
effectiveSize(herit)

# get the mean from the posterior disribution of heritability
mean(herit)

# get confidence interval for heritability
HPDinterval(herit)

#plot the trace of heritability, should not be any pattern
plot(herit)


#######################################################
# Model 1101: Trait = Ptype1, Batch fixed, Additive effects and maternal effects are random

#needed library for running genearlized linear mixed models
library(MCMCglmm)

#Clear memory if needed
rm(list=ls()) 

load("Toy4.rat")
load("Toy4ped.rat")

# set up parameters for the priors
# Fixed effects are automatically set by MCMCglmm to follow a normal distribution and do not need to be specified
# G is for specified random effects (additive, maternal, etc). 
# Here we set for weak priors that will be used in an inverse gamma disribution automatically set by MCMCglmm
# R is residual effects for each specified random effect and follows same rules as G
prior1101 <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002), G2 = list(V = 1, nu = 0.002)))

#Set up inverse of Numerator Relationship Matrix using MCMCglmm inverseA command
invA <- inverseA(Toy4ped) 

# model statement
model1101 <- MCMCglmm(Ptype1 ~ 1 + Batch,               #intercept (the 1) and Batch are the fixed effect
                     random = ~animal + Damid,          #additive (animal) and maternal (Damid) are random effects
                     family = "gaussian",               #phenotype has gaussian distribution
                     prior = prior1101,                 #call the priors parameters defined above
                     data = Toy4,                       #call the data file
                     nitt = 100000,                     #number of MCMC iterations 
                     burnin = 10000,                    #number of iterations for burnin
                     thin = 1000,                       #sampling interval
                     ginverse=list(animal=invA)$invA)   #call the inverse of the NRM
                     

#save model as an R object so we can access it later
save(model1101, file = "model1101.obj") #dont forget to change the path#

#load model file
load(file = "model1101.obj") #dont forget to change the path#

#plot trace and density of fixed effects; should be no trend in trace
plot(model1101$Sol)
#plot trace and density of random (additive, maternal and residual (=environmental)) variances; should be no trend in trace
plot(model1101$VCV)

#examine autocorrelation of fixed effects
autocorr.diag(model1101$Sol)
#examine autocorrelation of random (additive, maternal and residual) variances
autocorr.diag(model1101$VCV)

#check effective population size for fixed effects; should be gt 1000
effectiveSize(model1101$Sol)
#check effective population size for random effects (additve, maternal and residual variances); should be gt 1000
effectiveSize(model1101$VCV)

#test of convergence, p should be greater than 0.05 for good convergence
heidel.diag(model1101$VCV)

#estimates of additive and residual variances
posterior.mode(model1101$VCV)

#summary of model; make sure to check DIC score (smaller is better)
summary(model1101)

#estimate posterior distribution of the heritability (animal vairance divided by animal + residual variances)
herit <- model1101$VCV[, "animal"]/(model1101$VCV[, "animal"] + model1101$VCV[, "Damid"] + model1101$VCV[, "units"])

#effective sample size for heritability should be gt 1000
effectiveSize(herit)

# get the mean from the posterior disribution of heritability
mean(herit)

# get confidence interval for heritability
HPDinterval(herit)

#plot the trace of heritability, should not be any pattern
plot(herit)



#######################################################

#Class Exercises: 
#Estimate variance ratio for maternal effects.
#Rerun the above with maternal genetic effects.

#######################################################
# Model 1201: Traits = Ptype1 and Ptype2; Fixed = Batch, Random = Additive 

#Clear memory if needed
rm(list=ls()) 

#needed library for running genearlized linear mixed models
library(MCMCglmm)

# read in pedigree file; see details of pedigree data structure in model 101

load("Toy4ped.Rat")
# summarze pedigree file
summary(Toy4ped)

# read in data file; see details of data file structure in model 101

load("Toy4.Rat")
summary(Toy4)

# set up parameters for the priors
# Fixed effects are automatically set by MCMCglmm to follow a normal distribution and do not need to be specified
# G is for specified random effects (additive, maternal, etc). 
# Here we set for weak priors that will be used in an inverse gamma disribution automatically set by MCMCglmm
# R is residual effects for each specified random effect and follows same rules as G
prior1201 <- list(R=list(V=diag(2)*(0.002/1.002),nu=1.002),
                  G=list(G1=list(V=diag(2)*(0.002/1.002),nu=1.002)))

#Set up inverse of Numerator Relationship Matrix using MCMCglmm inverseA command
invA <- inverseA(Toy4ped) 

# model statement
model1201 <- MCMCglmm(cbind(Ptype1,Ptype2)~trait-1,                         #cbind combines the two traits into a matrix, intercept not fit, no fixed effects
                     random = ~us(trait):animal,                            #random effects, see below for explanation
                     rcov=~us(trait):units,                                 #residual effects, see below for explanation
                     family = c("gaussian","gaussian"),                     #both phenotypes have gaussian distribution
                     prior = prior1201,                                     #call the priors parameters
                     data = Toy4,                                           #call the data file
                     nitt = 100000,                                         #number of MCMC iterations: NORMALLY NEED ~10,000,000 for a good run and convergence. here is only 100,000 (1/100 optimal, for time's sake)
                     burnin = 10000,                                        #number of iterations for burnin
                     thin = 1000,                                           #sampling interval
                     ginverse=list(animal=invA)$invA)                       #call the inverse of the NMR

# A note on the random statement in a multivariate model
#
# The term random=~us(trait):animal means for the additive effect use the following structure for the G matrix:
#    V1    COV12      where V1 = variance of trait 1 and COV12 = covariance of traits 1 and 2
#    COV12 V2         where V2 = variance of trait 2 and COV12 = covariance of traits 1 and 2
# which is usually what we want.  

# If we want to set a G matrix with covariances of 0 we would put this term into our model
# random=~idh(trait):animal which means for the additive effect use the following structure for the G matrix
#    V1  0            where V1 and V2 are the variances of traits 1 and 2 respectively
#    0   V2           and the covariances are set to zero
#
# Finally we can put the following term in the random statement
# random=~animal
# which means use the following matrix
#   V  V
#   V  V
# i.e., all the variances and covariances are assumed to be equal
# The same logic and syntax apply to the rcov command which estimates the residual variances and covariances


#save model as an R object so we can access it later if needed
save(model1201, file = "model1201.obj")

#load model if needed
load("model1201.obj")

#plot trace and density of fixed effects; should be no trend in trace
plot(model1201$Sol)
#plot trace and density of random (additive and residual) variances; should be no trend in trace
plot(model1201$VCV)

#examine autocorrelation of fixed effects
autocorr.diag(model1201$Sol)
#examine autocorrelation of random (additive and residual) variances
autocorr.diag(model1201$VCV)

#check effective population size for fixed effects; should be gt 1000
effectiveSize(model1201$Sol)
#check effective population size for random effects (additve and residual variances); should be gt 1000
effectiveSize(model1201$VCV)

#test of convergence, p should be greater than 0.05 for good convergence
heidel.diag(model1201$VCV)

#estimates of additive and residual variances
posterior.mode(model1201$VCV)

#summary of model; make sure to check DIC score (smaller is better)
summary(model1201)

#estimate posterior distribution of the heritability for trait = Ptype1
heritPtype1 <- model1201$VCV[, "traitPtype1:traitPtype1.animal"]/(model1201$VCV[, "traitPtype1:traitPtype1.animal"] + model1201$VCV[, "traitPtype1:traitPtype1.units"])
#effective sample size for heritability of trait 1 (should be gt 1000)
effectiveSize(heritPtype1)
# get the mean from the posterior disribution of heritability 
mean(heritPtype1)
# get confidence interval for heritability
HPDinterval(heritPtype1)
#plot the trace of heritability, should not be any pattern
plot(heritPtype1)


#estimate posterior distribution of the heritability for trait = Ptype2
heritPtype2 <- model1201$VCV[, "traitPtype2:traitPtype2.animal"]/(model1201$VCV[, "traitPtype2:traitPtype2.animal"] + model1201$VCV[, "traitPtype2:traitPtype2.units"])
#effective sample size for heritability of Ptype2 (should be gt 1000)
effectiveSize(heritPtype2)
# get the mean from the posterior disribution of heritability
mean(heritPtype2)
# get confidence interval for heritability
HPDinterval(heritPtype2)
#plot the trace of heritability, should not be any pattern
plot(heritPtype2)


#estimate posterior distribution of the genetic correlation between Ptype1 and Ptype2 
GenCorrPtype1Ptype2 <-model1201$VCV[, "traitPtype1:traitPtype2.animal"]/sqrt(model1201$VCV[, "traitPtype1:traitPtype1.animal"]*model1201$VCV[, "traitPtype2:traitPtype2.animal"])
#effective sample size for genetic correlation (should be gt 1000)
effectiveSize(GenCorrPtype1Ptype2)
#mean of posterior distribution of genetic correlation
mean(GenCorrPtype1Ptype2)
#get confidence interval for genetic correlation
HPDinterval(GenCorrPtype1Ptype2)
#plot the trace of genetic correlation, should not be any pattern
plot(GenCorrPtype1Ptype2)

#######################################################
#Class exercises
#Run the above with Batch as a Fixed effect
#Run the above with maternal effects

#######################################################


#######################################################
# Model 2001: Traits = Ptype1, Ptype2, Ptype3; Fixed = Batch, Random = Additive 

#Clear memory if needed
rm(list=ls()) 

#needed library for running genearlized linear mixed models
library(MCMCglmm)

# read in pedigree file; see details of pedigree data structure in model 101

load("Toy4ped.Rat")
# summarze pedigree file
summary(Toy4ped)

# read in data file; see details of data file structure in model 101

load("Toy4.Rat")
summary(Toy4)

# set up parameters for the priors
# Fixed effects are automatically set by MCMCglmm to follow a normal distribution and do not need to be specified
# G is for specified random effects (additive, maternal, etc). 
# Here we set for weak priors that will be used in an inverse gamma disribution automatically set by MCMCglmm
# R is residual effects for each specified random effect and follows same rules as G
prior2001 <- list(R=list(V=diag(3)*(0.002/1.002),nu=1.002),
                  G=list(G1=list(V=diag(3)*(0.002/1.002),nu=1.002)))

#Set up inverse of Numerator Relationship Matrix using MCMCglmm inverseA command
invA <- inverseA(Toy4ped) 

# model statement
model2001 <- MCMCglmm(cbind(Ptype1,Ptype2,Ptype3)~trait-1 + Batch,           #cbind combines the three traits into a matrix, intercept not fit, Batch is fixed
                      random = ~us(trait):animal,                            #random effects, see below for explanation
                      rcov=~us(trait):units,                                 #residual effects, see below for explanation
                      family = c("gaussian","gaussian","gaussian"),          #all phenotypes have gaussian distribution
                      prior = prior2001,                                     #call the priors parameters
                      data = Toy4,                                           #call the data file
                      nitt = 100000,                                         #number of MCMC iterations
                      burnin = 10000,                                        #number of iterations for burnin
                      thin = 1000,                                           #sampling interval
                      ginverse=list(animal=invA)$invA)                       #call the inverse of the NMR

# For notes on the random statement in a multivariate model see Model 1201


#save model as an R object so we can access it later if needed
save(model2001, file = "model2001.obj")

#load model if needed
load("model2001.obj")

#plot trace and density of fixed effects; should be no trend in trace
plot(model2001$Sol)
#plot trace and density of random (additive and residual) variances; should be no trend in trace
plot(model2001$VCV)

#examine autocorrelation of fixed effects
autocorr.diag(model2001$Sol)
#examine autocorrelation of random (additive and residual) variances
autocorr.diag(model2001$VCV)

#check effective population size for fixed effects; should be gt 1000
effectiveSize(model2001$Sol)
#check effective population size for random effects (additve and residual variances); should be gt 1000
effectiveSize(model2001$VCV)

#test of convergence, p should be greater than 0.05 for good convergence
heidel.diag(model2001$VCV)

#estimates of additive and residual variances
posterior.mode(model2001$VCV)

#summary of model; make sure to check DIC score (smaller is better)
summary(model2001)

#estimate posterior distribution of the heritability for trait = Ptype1
heritPtype1 <- model2001$VCV[, "traitPtype1:traitPtype1.animal"]/(model2001$VCV[, "traitPtype1:traitPtype1.animal"] + model2001$VCV[, "traitPtype1:traitPtype1.units"])
#effective sample size for heritability of trait 1 (should be gt 1000)
effectiveSize(heritPtype1)
# get the mean from the posterior disribution of heritability 
mean(heritPtype1)
# get confidence interval for heritability
HPDinterval(heritPtype1)
#plot the trace of heritability, should not be any pattern
plot(heritPtype1)


#estimate posterior distribution of the heritability for trait = Ptype2
heritPtype2 <- model2001$VCV[, "traitPtype2:traitPtype2.animal"]/(model2001$VCV[, "traitPtype2:traitPtype2.animal"] + model2001$VCV[, "traitPtype2:traitPtype2.units"])
#effective sample size for heritability of PC2 (should be gt 1000)
effectiveSize(heritPtype2)
# get the mean from the posterior disribution of heritability
mean(heritPtype2)
# get confidence interval for heritability
HPDinterval(heritPtype2)
#plot the trace of heritability, should not be any pattern
plot(heritPtype2)


#estimate posterior distribution of the heritability for trait = Ptype3
heritPtype3 <- model2001$VCV[, "traitPtype3:traitPtype3.animal"]/(model2001$VCV[, "traitPtype3:traitPtype3.animal"] + model2001$VCV[, "traitPtype3:traitPtype3.units"])
#effective sample size for heritability of trait 1 (should be gt 1000)
effectiveSize(heritPtype3)
# get the mean from the posterior disribution of heritability 
mean(heritPtype3)
# get confidence interval for heritability
HPDinterval(heritPtype3)
#plot the trace of heritability, should not be any pattern
plot(heritPtype3)


#estimate posterior distribution of the genetic correlation between Ptype1 Ptype2
GenCorrPtype1Ptype2 <-model2001$VCV[, "traitPtype1:traitPtype2.animal"]/sqrt(model2001$VCV[, "traitPtype1:traitPtype1.animal"]*model2001$VCV[, "traitPtype2:traitPtype2.animal"])
#effective sample size for genetic correlation (should be gt 1000)
effectiveSize(GenCorrPtype1Ptype2)
#mean of posterior distribution of genetic correlation
mean(GenCorrPtype1Ptype2)
#get confidence interval for genetic correlation
HPDinterval(GenCorrPtype1Ptype2)
#plot the trace of genetic correlation, should not be any pattern
plot(GenCorrPtype1Ptype2)


#estimate posterior distribution of the genetic correlation between Ptype1 and Ptype3 
GenCorrPtype1Ptype3 <-model2001$VCV[, "traitPtype1:traitPtype3.animal"]/sqrt(model2001$VCV[, "traitPtype1:traitPtype1.animal"]*model2001$VCV[, "traitPtype3:traitPtype3.animal"])
#effective sample size for genetic correlation (should be gt 1000)
effectiveSize(GenCorrPtype1Ptype3)
#mean of posterior distribution of genetic correlation
mean(GenCorrPtype1Ptype3)
#get confidence interval for genetic correlation
HPDinterval(GenCorrPtype1Ptype3)
#plot the trace of genetic correlation, should not be any pattern
plot(GenCorrPtype1Ptype3)


#estimate posterior distribution of the genetic correlation between Ptype2 and Ptype3 
GenCorrPtype2Ptype3 <-model2001$VCV[, "traitPtype2:traitPtype3.animal"]/sqrt(model2001$VCV[, "traitPtype2:traitPtype2.animal"]*model2001$VCV[, "traitPtype3:traitPtype3.animal"])
#effective sample size for genetic correlation (should be gt 1000)
effectiveSize(GenCorrPtype2Ptype3)
#mean of posterior distribution of genetic correlation
mean(GenCorrPtype2Ptype3)
#get confidence interval for genetic correlation
HPDinterval(GenCorrPtype2Ptype3)
#plot the trace of genetic correlation, should not be any pattern
plot(GenCorrPtype2Ptype3)




