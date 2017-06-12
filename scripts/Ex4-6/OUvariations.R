###################################################
### chunk number 1: 
###################################################
require(ouch)

regimes <- read.csv('data/Ex4-6/regimes.csv',row.names=1)    # regimes
regimes$BM <- "BM"
ssd <- read.csv('data/Ex4-6/ssd.data.csv',row.names=1)    # tree + body size data
otree <- with(ssd,ouchtree(nodes,ancestors,times,labels)) #ouch tree
plot(otree)
xdata <- log(ssd[c('fSVL','mSVL')])
names(xdata) <- paste('log',names(xdata),sep='.')
nreg <- length(regimes)


###################################################
### chunk number 2: 
###################################################
alpha.guess <- c(1, 0, 1)
sigma.guess <- c(1, 1, 1)


###################################################
### chunk number 3: 
###################################################
tic <- Sys.time()
h.OU6TW <- hansen(
       data=xdata,
       tree=otree,
       regimes=regimes["OU6TW"],
       sqrt.alpha =alpha.guess,
       sigma=sigma.guess,
       method="Nelder-Mead",
       maxit=3000,
       reltol=1e-12
       )
toc <- Sys.time()
print(toc-tic)
summary(h.OU6TW)


###################################################
### chunk number 4: 
###################################################

#increased max interations and relative tolerance to improve chances of convergence
tic <- Sys.time()
h.OU7LP <- hansen(
       data=xdata,
       tree=otree,
       regimes=regimes["OU7LP"],
       sqrt.alpha =alpha.guess,
       sigma=sigma.guess,
       method="subplex",      # subplex is a much more robust optimizer
       maxit=20000,
       reltol=1e-12
       )
toc <- Sys.time()
print(toc-tic)
summary(h.OU7LP)


###################################################
### chunk number 5: 
###################################################
brown.fit <- brown(
               data=xdata,
               tree=otree
             )
summary(brown.fit)             

plot(otree, regimes=regimes[c("BM", "OU6TW", "OU7LP")])

h <- list(brown.fit, h.OU6TW, h.OU7LP)
names(h) <- c("BM", "OU6TW", "OU7LP")
h.ic <- sapply( h, function(x) unlist( summary(x)[c("aic", "aic.c", "sic", "dof")]))
print( h.ic, digits = 3)

