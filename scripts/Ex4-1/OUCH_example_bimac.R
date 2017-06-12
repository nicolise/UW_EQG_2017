require(ouch)
data(bimac)
bimac

#a builtin data set
rownames(bimac) <- bimac$node

tree <- with( bimac, ouchtree(node, ancestor, time=time/max(time), species))

plot(tree)
str(tree) #not that pretty. S4 data objects, can't change them directly

#coloring tree based on the different phenotypic state hypotheses
#this includes all nodes, even the most basal nodes (less information on these, of course) -- but you *could* add in fossil data!
#prior to this, must build a tree with hypothesized time from root to each node
#tips are a special type of node

#everything matches based on node name
plot(tree, regimes=bimac[c("OU.1", "OU.3", "OU.4", "OU.LP")], lwd=3)

#input quantitative size measurements!
#log transform: animals tend to grow nonlinearly (bigger, grow faster)
size <- log(bimac$size) 
names(size) <- bimac$node

#annotate your hypotheses
h1 <- brown( size, tree) # BM
h2 <- hansen( size, tree, regimes=bimac['OU.1'], sqrt.alpha=1, sigma=1)
h3 <- hansen( size, tree, regimes=bimac['OU.3'], sqrt.alpha=1, sigma=1)
h4 <- hansen( size, tree, regimes=bimac['OU.4'], sqrt.alpha=1, sigma=1)
h5 <- hansen( size, tree, regimes=bimac['OU.LP'], sqrt.alpha=1, sigma=1)

coef(h5)
#grab just a subset of values from the summary -- into a vector
unlist( summary(h5)[c("aic", "aic.c", "sic", "dof")])

h <- list(h1, h2, h3, h4, h5)
names(h) <- c("BM", "OU.1", "OU.3", "OU.4", "OU.LP")

#unlist: lose dataframe properties
#sapply: apply function once per element in vector, turn it into a matrix
h.ic <- sapply( h, function(x) unlist( summary(x)[c("aic", "aic.c", "sic", "dof")]))
print( h.ic, digits = 3)

#a parametric bootstrap -- 10 sets of simulated data
h5.sim <- simulate(object=h5, nsim=10)
#update your best model
summary(update(object = h5, data= h5.sim[[1]]))

h5.sim.fit <- lapply(h5.sim, function(x) update(h5, x))
bootstrap (object = h5, nboot=10)

#want to bootstrap model selection too: want to bootstrap your best model, and then simulate the alternative hypotheses
plot(tree, node.names=T)
ou.lp <- paint(tree, subtree=c("1"="medium","9"="large","2"="small") )
plot(tree, regimes=ou.lp, node.names=T)

ou.lp <- paint(tree, subtree=c("1"="medium","9"="large","2"="small"), branch=c("38"="large","2"="medium"))

plot(tree, regimes=ou.lp, node.names=T)

ou.clades <- paint( tree, subtree=c("1"="A","7"="B", "8"="C"), branch=c("8"="C", "7"="C", "1"="A"))
plot(tree, regimes=ou.clades, node.names=T)

h6 <- hansen(log(bimac['size']),tree, regimes=ou.clades,sqrt.alpha=1,sigma=1)
h <- append(h, h6)         # append (add on) new model results to our list h
names(h)[length(h)] <- "OU.clades"    # add the name of the new model
names(h)
h.ic <- sapply( h, function(x) unlist(summary(x)[c('aic', 'aic.c', 'sic', 'dof')]) ) 
print( h.ic, digits = 3)
