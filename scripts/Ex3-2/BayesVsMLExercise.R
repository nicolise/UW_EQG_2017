#Bayes vs ML exercise
#Brian O'Meara, 6 June 2017


#let's look at a log likelihood plot
xvals<-seq(from=-1, to=1, length.out=1000)
plot(xvals, log(sapply(xvals, dnorm, mean=0, sd=4)), type="l", bty="n", xlab="value", ylab="Log likelihood",ylim=c(-2,2))
abline(h=0, lty="dotted")
sd.vector <- seq(from=1, to=0.05, length.out=5)
#dnorm = likelihood of those data under a normal distribution
#iterates over the vector. For each element, draws a normal curve with mean=0 and sd=value under that element.
for (iterator in sequence(length(sd.vector))) {
  lines(xvals, log(sapply(xvals, dnorm, mean=0, sd=sd.vector[iterator])))
}

#What? A positive log likelihood? Let's look at that in non-log space:
plot(xvals, sapply(xvals, dnorm, mean=0, sd=0.05), type="l", bty="n", xlab="value", ylab="Likelihood")
#How can we have a likelihood greater than 1?


#set up for next visualization
par(mfcol=c(3,3))

#Let's go back to coin flipping in R. Let's have three hypotheses:
#p(D)
prob.heads.hypotheses <- c(0.1, 0.5, 0.8)

#
#p(H)
prior.on.hypotheses <- c(0.05, 0.9, 0.05) #add your own here. What is the prob of c(0.1, 0.5, etc.
#this makes sure priors sum to 1, in case you are bad at addition
prior.on.hypotheses <- prior.on.hypotheses/sum(prior.on.hypotheses) #priors should sum to 1

#next assume our next flips we get 2/3 heads
num.flips <- 3
num.heads <- round(num.flips*2/3)
num.heads <- round(num.flips*1/2)

#how to calculate the likelihood?

#hint: ?dbinom
#dbinom(x=data, size, prob)
dbinom(num.heads, size=num.flips, prob=prob.heads.hypotheses[1])

#p(D|H)
likelihood.data.given.hypotheses <- c(dbinom(num.heads, size=num.flips, prob=prob.heads.hypotheses[1]), dbinom(num.heads, size=num.flips, prob=prob.heads.hypotheses[2]),dbinom(num.heads, size=num.flips, prob=prob.heads.hypotheses[3]))

#Now, write a function for Bayes Rule. Assume that the three hypotheses are the only ones possible

#store the posteriors as the vector posterior.on.hypotheses
#posterior.on.hypotheses<- likelihood.data.given.hypotheses*prior.on.hypotheses/prob.heads.hypotheses ... not quite

prob.data.given.hypothesis<-likelihood.data.given.hypotheses*prior.on.hypotheses
posterior.on.hypotheses <- prob.data.given.hypothesis/sum(prob.data.given.hypothesis)
  

barplot(prior.on.hypotheses, main="Prior")
barplot(likelihood.data.given.hypotheses, main="Likelihood")
barplot(posterior.on.hypotheses, main="Posterior")

#Now, try again with different priors
#And try again with different sample size (num.heads)

#Let's model evolution with a trend.
#We're going to go very basic for this. Our tree is a polytomy.
ntaxa<-10
sd.along.branch <- 0.3
trend.direction<-rnorm(1, 1, sd=0.5)
starting.state<-rnorm(1, sd=0.1)
tips<-rnorm(20,mean=trend.direction+starting.state, sd=sd.along.branch)
par(mfcol=c(1,1))
plot(x=c(0,1), y=range(c(starting.state, tips)), bty="n", xlab="time", ylab="state", type="n")
for (i in sequence(ntaxa)) {
  lines(c(0,1), c(starting.state, tips[i]))
}

#Now:
#Step 1: Program way to calculate likelihood of the observed tip.data.
GetLogLikelihood <- function(starting.state, trend.direction, sd.along.branch, tips) {
  individual.log.likelihoods<- ----SOMETHING----- #using log
  overall.log.likelihood <- -----SOMETHING------
  return(overall.log.likelihood)
}

#Does anyone see a problem at this point?


if(!require(akima)) {
	install.packages("akima", repos="http://cran.us.r-project.org")
}
library(akima)

nreps=400
starting.vector <- rnorm(400, starting.state, .2)
trend.vector <- rnorm(400, trend.direction, .2)
likelihood.values<-rep(NA,nreps)
for (i in sequence(nreps)) {
  likelihood.values[i] <- GetLogLikelihood(starting.vector[i], trend.vector[i], sd.along.branch, tips)
}

plot(starting.vector + trend.vector, likelihood.values, xlab="Trend + Starting", ylab="Log Likelihood", bty="n", pch=".")

#looks good, right?

rescaled.likelihood<-max(likelihood.values) - likelihood.values


interpolated.points<-interp(x=starting.vector, y=trend.vector, z= rescaled.likelihood, linear=FALSE, extrap=TRUE, xo=seq(min(starting.vector), max(starting.vector), length = 400), yo=seq(min(trend.vector), max(trend.vector), length = 400))

contour(interpolated.points, xlim=range(starting.vector),ylim=range(trend.vector), xlab="Starting state", ylab="Trend", levels=c(2,5,10),add=FALSE,lwd=1, bty="n", asp=1)

#How would a prior affect this?


#AIC vs AICc vs BIC

# First, make a function to go from log likelihood, number of parameters, and number of datapoints to AIC, AICc, and BIC

GetICScores <- function(lnl, n, k, not.just.penalty=FALSE) {
  AIC = (-2 * lnl)*not.just.penalty + 2 * k
  AICc = (-2 * lnl)*not.just.penalty + 2 * k * n / (n - k - 1)
  BIC = (-2 * lnl)*not.just.penalty + log(n) * k
  result.vector <- c(AIC, AICc, BIC)
  names(result.vector) <- c("AIC", "AICc", "BIC")
  return(result.vector)
}

# Now generate 1000 data points under a gamma distribution
sim.points <- rgamma(n=1000, shape=3, rate=0.1)

# Always, always, always look at your data:

plot(density(sim.points))

print(quantile(sim.points, prob=seq(from=0, to=1, length.out=11)))

# Now calculate the likelihoods for each, under a normal distribution and an exponential distribution. Why use a distribution different from the one we used? b/c in life we often use a different distribution from the true one. For simplicity, we're going to do rough estimates of the values, not truly optimize.

sim.mean <- mean(sim.points)
sim.sd <- sd(sim.points)

# these give the log likelihoods for each point. We can sum them to get the log likelihood for all together
lnl.normal <- dnorm(sim.points, mean=sim.mean, sd=sim.sd, log=TRUE) #uses two params. Normally you'd optimize them; here, we're just plugging in values that are pretty good
lnl.exp<- dexp(sim.points, rate=1/sim.mean, log=TRUE) #uses one param; see comment above about optimizing

GetPlottingPoints <- function(n.to.sample, lnl.normal, lnl.exp, k.normal=2, k.exp=1, not.just.penalty=TRUE) {
   normal.ic <- GetICScores(lnl = n.to.sample*mean(lnl.normal), n=n.to.sample, k=k.normal, not.just.penalty = not.just.penalty)
   exp.ic <-  GetICScores(lnl = n.to.sample*mean(lnl.exp), n=n.to.sample, k=k.exp, not.just.penalty = not.just.penalty)
   results <- c(n.to.sample*mean(lnl.normal), normal.ic, n.to.sample*mean(lnl.exp), exp.ic)
   names(results) <- c("lnl_normal", "AIC_normal", "AICc_normal", "BIC_normal", "lnl_exponential", "AIC_exponential", "AICc_exponential", "BIC_exponential")
   return(results)
}

point.counts <- seq(from=4, to=20, by=1)
scores.penalties <- sapply(point.counts, GetPlottingPoints, lnl.normal=lnl.normal, lnl.exp=lnl.exp, not.just.penalty=FALSE)
scores.full <- sapply(point.counts, GetPlottingPoints, lnl.normal=lnl.normal, lnl.exp=lnl.exp, not.just.penalty=TRUE)

colnames(scores.penalties) <- point.counts
colnames(scores.full) <- point.counts

par(mfcol=c(3,3))
method.names <- c("AIC", "AICc", "BIC")
models <- c("normal", "exponential")
cols <- c("red", "black")
for (method.type in sequence(3)) {
  plot(x=range(as.numeric(colnames(scores.full))), y=range(scores.full[!grepl("lnl", rownames(scores.full)),]), xlab="Number of points", ylab="Overall IC score", bty="n", type="n", main=method.names[method.type])
  relevant.rows <- grepl(paste0(method.names[method.type],"_"), rownames(scores.full))
  scores.penalties.local <- scores.penalties[relevant.rows,]
  scores.full.local <- scores.full[relevant.rows,]
  for (model.type in sequence(length(models))) {
    lines(x=as.numeric(colnames(scores.full)), y=scores.full.local[grepl(models[model.type],rownames(scores.penalties.local)),], col=cols[model.type])
  }
}

for (method.type in sequence(3)) {
  plot(x=range(as.numeric(colnames(scores.full))), y=c(-1,1)*3, xlab="Number of points", ylab="∆score: When negative, normal is better", bty="n", type="n", main=method.names[method.type])
  relevant.rows <- grepl(paste0(method.names[method.type],"_"), rownames(scores.full))
  scores.full.local <- scores.full[relevant.rows,]
  exp.scores <- scores.full.local[grepl("exponential",rownames(scores.full.local)),]
  norm.scores <- scores.full.local[grepl("normal",rownames(scores.full.local)),]
  lines(x=as.numeric(colnames(scores.full)), y=norm.scores - exp.scores, col=cols[model.type], lwd=2)
  abline(h=0, col="purple")

}


for (method.type in sequence(3)) {
  plot(x=range(as.numeric(colnames(scores.full))), y=c(0,max(scores.penalties[!grepl("lnl", rownames(scores.penalties)),])), xlab="Number of points", ylab="Penalty term", bty="n", type="n", main=method.names[method.type])
  relevant.rows <- grepl(paste0(method.names[method.type],"_"), rownames(scores.full))
  scores.penalties.local <- scores.penalties[relevant.rows,]
  scores.full.local <- scores.full[relevant.rows,]
  for (model.type in sequence(length(models))) {
    lines(x=as.numeric(colnames(scores.full)), y=scores.penalties.local[grepl(models[model.type],rownames(scores.penalties.local)),], col=cols[model.type], lwd=2)
  }
}

# What does this mean? When is exponential better than normal -- for AIC? AICc? BIC?
