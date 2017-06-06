
# This R script was written by Patrick Carter in August 2016 for the NIMBioS Evolutionary Genetics course
# The purpose of this script is to highlight code needed to run quantitative genetic analyses in R using MCMCglmm
# The data to be analyzed are body size traits for farmed rainbow trout: head length, body length, body depth, body mass
# The purpose of the project from which these data were drawn was to predict the response to selection on body depth
# and to predict correlated responses to this selection on the other phenotypes

#IMPORTANT NOTE: You will need to change the path for opening data files to whatever is appropriate to your computer


#######################################################
# Model 101: Trait = Head Length, Batch fixed, Additive effects only random effect

#Clear memory if needed
rm(list=ls()) 

#needed library for running genearlized linear mixed models
library(MCMCglmm)

# read in pedigree file 

# Important information about the structure and construction of pedigree files to be used by MCMCglmm
# pedigree files contain identification information for measured individuals and their mothers (Dam) and fathers (Sire) and other ancestors
# this file may contain many generations of individuals for whom phenotypes were never measured
# note that offspring always must have larger id numbers than their parents, parents must have larger id numbers than grandparents, etc

# make sure missing values for Sire and Dam identification numbers are coded as NA
# make sure that dam and sire id numbers are coded as factors
# the identifiation variable for the individuals with phenotypic information in the data set MUST called animal and it should remain a numeric variable
# the pedigree file should be sorted by the variable animal

# individuals may appear in the pedgiree as animal
# and they may also appear as a Sire or Dam

# this commmand from the MasterBayes library completes a pedigree with missing information for some sires and dams 
# by adding the generation in which all Dams and Sires were unknown; you will need to do this when making your pedigree file:
# MyPedigree<-insertPed(MyPedigree, founders=NULL)

FishPed<-read.delim ("C:\\Users\\Patrick Carter\\Documents\\Teaching\\Nimbios16\\FishPed.txt")  #dont forget to change the path#
#Make sure Dam and Sire are factors
FishPed$Dam <- as.factor(FishPed$Dam)
FishPed$Sire <- as.factor(FishPed$Sire)
# summarze pedigree file
summary(FishPed)

# read in data file

# make sure missing values for Sire and Dam id numbers are coded as NA
# make sure that dam and sire id numbers are coded as factors
# the identifiation variable for individuals with phenotypic data MUST called animal and it should remain a numeric variable
# drop cases of individuals without phenotypic data 
# file should be sorted by the variable animal

FishHBDM<-read.delim ("C:\\Users\\Patrick Carter\\Documents\\Teaching\\Nimbios16\\FishHBDM.txt")  #dont forget to change the path#
#Make sure Dam and Sire are factors
FishHBDM$Dam <- as.factor(FishHBDM$Dam)
FishHBDM$Sire <- as.factor(FishHBDM$Sire)
#summarize data file (Note: before this point you should already have thoroughly graphed and examined your data)
summary(FishHBDM)

#check the phenotypic variance/covariance matrix for the 4 phenotypes (head length, body length, body depth, body mass)
#you should also have graphically inspected your data before you run MCMCglmm
TempData <- subset(FishHBDM, select = c(Head,Body,Depth))
VarCovMatrix = cov(TempData)

# set up parameters for the priors
# Fixed effects are automatically set by MCMCglmm to follow a normal distribution and do not need to be specified
# G is for specified random effects (additive, maternal, etc). 
# Here we set for weak priors that will be used in an inverse gamma disribution automatically set by MCMCglmm
# R is residual effects for each specified random effect and follows same rules as G
prior101 <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002)))

#Set up inverse of Numerator Relationship Matrix using MCMCglmm inverseA command
invA <- inverseA(FishPed) 

# model statement
model101 <- MCMCglmm(Head ~ 1 + Batch,                  #intercept (the 1) and Batch are the fixed effect
                     random = ~animal,                  #additive effects (animal) the only random effect
                     family = "gaussian",               #phenotype has gaussian distribution
                     prior = prior101,                  #call the priors parameters defined above
                     data = FishHBDM,                   #call the data file
                     nitt = 1000000,                    #number of MCMC iterations 
                     burnin = 10000,                    #number of iterations for burnin
                     thin = 1000,                       #sampling interval
                     ginverse=list(animal=invA)$invA)   #call the inverse of the NRM

#save model as an R object so we can access it later
save(model101, file = "C:\\Users\\Patrick Carter\\Documents\\Teaching\\Nimbios16\\model101.obj") #dont forget to change the path#

#load model file
load(file = "C:\\Users\\Patrick Carter\\Documents\\Teaching\\Nimbios16\\model101.obj") #dont forget to change the path#

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

#Class Exercises: 
#Run the above model without batch as a fixed effect. 
#Run the above with maternal effects as a another random effect. 

#######################################################
# Model 1201: Traits = Head Length and Body Length; Fixed = Batch, Random = Additive 

#Clear memory if needed
rm(list=ls()) 

#needed library for running genearlized linear mixed models
library(MCMCglmm)

# read in pedigree file; see details of pedigree data structure in model 101

load("C:\\Users\\Patrick Carter\\Documents\\Teaching\\Nimbios16\\FishPed.Rat")
# summarze pedigree file
summary(FishPed)

# read in data file; see details of data file structure in model 101

load("C:\\Users\\Patrick Carter\\Documents\\Teaching\\Nimbios16\\FishHBDM.Rat")
summary(FishHBDM)

#check the phenotypic variance/covariance matrix for the 4 phenotypes (head length, body length, body depth, body mass)
TempData <- subset(FishHBDM, select = c(Head,Body,Depth))
VarCovMatrix = cov(TempData)

# set up parameters for the priors
# Fixed effects are automatically set by MCMCglmm to follow a normal distribution and do not need to be specified
# G is for specified random effects (additive, maternal, etc). 
# Here we set for weak priors that will be used in an inverse gamma disribution automatically set by MCMCglmm
# R is residual effects for each specified random effect and follows same rules as G
prior1201 <- list(R=list(V=diag(2)*(0.002/1.002),nu=1.002),
                  G=list(G1=list(V=diag(2)*(0.002/1.002),nu=1.002)))

#Set up inverse of Numerator Relationship Matrix using MCMCglmm inverseA command
invA <- inverseA(FishPed) 

# model statement
model1201 <- MCMCglmm(cbind(Head,Body)~trait-1,                             #cbind combines the two traits into a matrix, intercept not fit, no fixed effects
                     random = ~us(trait):animal,                            #random effects, see below for explanation
                     rcov=~us(trait):units,                                 #residual effects, see below for explanation
                     family = c("gaussian","gaussian"),                     #both phenotypes have gaussian distribution
                     prior = prior1201,                                     #call the priors parameters
                     data = FishHBDM,                                       #call the data file
                     nitt = 100000,                                         #number of MCMC iterations
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
save(model1201, file = "C:\\Users\\Patrick Carter\\Documents\\Teaching\\Nimbios16\\model1201.obj")

#load model if needed
load("C:\\Users\\Patrick Carter\\Documents\\Teaching\\Nimbios16\\model1201.obj")

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

#estimate posterior distribution of the heritability for trait = Head Length
heritHead <- model1201$VCV[, "traitHead:traitHead.animal"]/(model1201$VCV[, "traitHead:traitHead.animal"] + model1201$VCV[, "traitHead:traitHead.units"])
#effective sample size for heritability of trait 1 (should be gt 1000)
effectiveSize(heritHead)
# get the mean from the posterior disribution of heritability 
mean(heritHead)
# get confidence interval for heritability
HPDinterval(heritHead)
#plot the trace of heritability, should not be any pattern
plot(heritHead)


#estimate posterior distribution of the heritability for trait = Body Length 
heritBody <- model1201$VCV[, "traitBody:traitBody.animal"]/(model1201$VCV[, "traitBody:traitBody.animal"] + model1201$VCV[, "traitBody:traitBody.units"])
#effective sample size for heritability of PC2 (should be gt 1000)
effectiveSize(heritBody)
# get the mean from the posterior disribution of heritability
mean(heritBody)
# get confidence interval for heritability
HPDinterval(heritBody)
#plot the trace of heritability, should not be any pattern
plot(heritBody)


#estimate posterior distribution of the genetic correlation between Head Length and Body Length 
GenCorrHeadBody <-model1201$VCV[, "traitHead:traitBody.animal"]/sqrt(model1201$VCV[, "traitHead:traitHead.animal"]*model1201$VCV[, "traitBody:traitBody.animal"])
#effective sample size for genetic correlation (should be gt 1000)
effectiveSize(GenCorrHeadBody)
#mean of posterior distribution of genetic correlation
mean(GenCorrHeadBody)
#get confidence interval for genetic correlation
HPDinterval(GenCorrHeadBody)
#plot the trace of genetic correlation, should not be any pattern
plot(GenCorrHeadBody)

#######################################################
#Class exercises
#Run the above with Batch as a Fixed effect
#Run the above with maternal effects

#######################################################


#######################################################
# Model 2001: Traits = Head Length and Body Length and Body Depth; Fixed = Batch, Random = Additive 

#Clear memory if needed
rm(list=ls()) 

#needed library for running genearlized linear mixed models
library(MCMCglmm)

# read in pedigree file; see details of pedigree data structure in model 101

load("C:\\Users\\Patrick Carter\\Documents\\Teaching\\Nimbios16\\FishPed.Rat")
# summarze pedigree file
summary(FishPed)

# read in data file; see details of data file structure in model 101

load("C:\\Users\\Patrick Carter\\Documents\\Teaching\\Nimbios16\\FishHBDM.Rat")
summary(FishHBDM)

#check the phenotypic variance/covariance matrix for the 4 phenotypes (head length, body length, body depth, body mass)
TempData <- subset(FishHBDM, select = c(Head,Body,Depth))
VarCovMatrix = cov(TempData)

# set up parameters for the priors
# Fixed effects are automatically set by MCMCglmm to follow a normal distribution and do not need to be specified
# G is for specified random effects (additive, maternal, etc). 
# Here we set for weak priors that will be used in an inverse gamma disribution automatically set by MCMCglmm
# R is residual effects for each specified random effect and follows same rules as G
prior2001 <- list(R=list(V=diag(3)*(0.002/1.002),nu=1.002),
                  G=list(G1=list(V=diag(3)*(0.002/1.002),nu=1.002)))

#Set up inverse of Numerator Relationship Matrix using MCMCglmm inverseA command
invA <- inverseA(FishPed) 

# model statement
model2001 <- MCMCglmm(cbind(Head,Body,Depth)~trait-1 + Batch,                #cbind combines the three traits into a matrix, intercept not fit, Batch is fixed
                      random = ~us(trait):animal,                            #random effects, see below for explanation
                      rcov=~us(trait):units,                                 #residual effects, see below for explanation
                      family = c("gaussian","gaussian","gaussian"),          #all phenotypes have gaussian distribution
                      prior = prior2001,                                     #call the priors parameters
                      data = FishHBDM,                                       #call the data file
                      nitt = 100000,                                         #number of MCMC iterations
                      burnin = 10000,                                        #number of iterations for burnin
                      thin = 1000,                                           #sampling interval
                      ginverse=list(animal=invA)$invA)                       #call the inverse of the NMR

# For notes on the random statement in a multivariate model see Model 1201


#save model as an R object so we can access it later if needed
save(model2001, file = "C:\\Users\\Patrick Carter\\Documents\\Teaching\\Nimbios16\\model2001.obj")

#load model if needed
load("C:\\Users\\Patrick Carter\\Documents\\Teaching\\Nimbios16\\model2001.obj")

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

#estimate posterior distribution of the heritability for trait = Head Length
heritHead <- model2001$VCV[, "traitHead:traitHead.animal"]/(model2001$VCV[, "traitHead:traitHead.animal"] + model2001$VCV[, "traitHead:traitHead.units"])
#effective sample size for heritability of trait 1 (should be gt 1000)
effectiveSize(heritHead)
# get the mean from the posterior disribution of heritability 
mean(heritHead)
# get confidence interval for heritability
HPDinterval(heritHead)
#plot the trace of heritability, should not be any pattern
plot(heritHead)


#estimate posterior distribution of the heritability for trait = Body Length 
heritBody <- model2001$VCV[, "traitBody:traitBody.animal"]/(model2001$VCV[, "traitBody:traitBody.animal"] + model2001$VCV[, "traitBody:traitBody.units"])
#effective sample size for heritability of PC2 (should be gt 1000)
effectiveSize(heritBody)
# get the mean from the posterior disribution of heritability
mean(heritBody)
# get confidence interval for heritability
HPDinterval(heritBody)
#plot the trace of heritability, should not be any pattern
plot(heritBody)


#estimate posterior distribution of the heritability for trait = Body Depth
heritDepth <- model2001$VCV[, "traitDepth:traitDepth.animal"]/(model2001$VCV[, "traitDepth:traitDepth.animal"] + model2001$VCV[, "traitDepth:traitDepth.units"])
#effective sample size for heritability of trait 1 (should be gt 1000)
effectiveSize(heritDepth)
# get the mean from the posterior disribution of heritability 
mean(heritDepth)
# get confidence interval for heritability
HPDinterval(heritDepth)
#plot the trace of heritability, should not be any pattern
plot(heritDepth)


#estimate posterior distribution of the genetic correlation between Head Length and Body Length 
GenCorrHeadBody <-model2001$VCV[, "traitHead:traitBody.animal"]/sqrt(model2001$VCV[, "traitHead:traitHead.animal"]*model2001$VCV[, "traitBody:traitBody.animal"])
#effective sample size for genetic correlation (should be gt 1000)
effectiveSize(GenCorrHeadBody)
#mean of posterior distribution of genetic correlation
mean(GenCorrHeadBody)
#get confidence interval for genetic correlation
HPDinterval(GenCorrHeadBody)
#plot the trace of genetic correlation, should not be any pattern
plot(GenCorrHeadBody)


#estimate posterior distribution of the genetic correlation between Head Length and Body Depth 
GenCorrHeadDepth <-model2001$VCV[, "traitHead:traitDepth.animal"]/sqrt(model2001$VCV[, "traitHead:traitHead.animal"]*model2001$VCV[, "traitDepth:traitDepth.animal"])
#effective sample size for genetic correlation (should be gt 1000)
effectiveSize(GenCorrHeadDepth)
#mean of posterior distribution of genetic correlation
mean(GenCorrHeadDepth)
#get confidence interval for genetic correlation
HPDinterval(GenCorrHeadDepth)
#plot the trace of genetic correlation, should not be any pattern
plot(GenCorrHeadDepth)


#estimate posterior distribution of the genetic correlation between Body Length and Body Depth 
GenCorrBodyDepth <-model2001$VCV[, "traitBody:traitDepth.animal"]/sqrt(model2001$VCV[, "traitBody:traitBody.animal"]*model2001$VCV[, "traitDepth:traitDepth.animal"])
#effective sample size for genetic correlation (should be gt 1000)
effectiveSize(GenCorrBodyDepth)
#mean of posterior distribution of genetic correlation
mean(GenCorrBodyDepth)
#get confidence interval for genetic correlation
HPDinterval(GenCorrBodyDepth)
#plot the trace of genetic correlation, should not be any pattern
plot(GenCorrBodyDepth)


#######################################################
#Class exercises
#Run the above with all 4 phenotypes

#######################################################
