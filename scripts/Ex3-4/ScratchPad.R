rm(list=ls())
nsteps=100
#random normal number generater
#default mean = 0, sd = 1
devs <- rnorm(100)

#for now, x is a dummy. 100 numbers in order
x <- c(0:100)
cbind(x,devs) 

x[1]
#how to get x2?
x[2] <- x[1] + devs[1]

#now generalize that
for (i in 1:nsteps){
  x[i+1] <- x[i] + devs[i]
}

myx <- as.data.frame(cbind(x,devs))
plot(myx$x, myx$devs)

#actual plot
plot(1:length(devs), devs, type = "n", col="red", 
     ylim = c(-max(devs)*30, max(devs)*30), 
     xlab="Time", ylab="Value", main="BM Simulation")
#simple version
x <- c(0:100)
devs <- rnorm(nsteps)
for (i in 1:nsteps) 
{
  x[i+1] <- x[i] + devs[i]
  lines(i:(i+1), x[i:(i+1)], col="red")
}

#now add sigma
#actual plot
plot(1:length(devs), devs, type = "n", col="red", 
     ylim = c(-max(devs)*30, max(devs)*30), 
     xlab="Time", ylab="Value", main="BM Simulation")
#simple version
x <- c(0:100)
devs <- rnorm(nsteps); sigma = 3
for (i in 1:nsteps) 
{
  x[i+1] <- x[i] + sigma*devs[i]
  lines(i:(i+1), x[i:(i+1)], col="red")
}

#make it a function!
#increase sigma: a drunker walk (bigger steps)
bmsim <- function(sigma=1){
  plot(1:length(devs), devs, type = "n", col="red", 
  ylim = c(-max(devs)*30, max(devs)*30), 
  xlab="Time", ylab="Value", main="BM Simulation")
  #simple version
  x <- c(0:100)
  devs <- rnorm(nsteps)
  for (i in 1:nsteps) 
  {
    x[i+1] <- x[i] + sigma*devs[i]
    lines(i:(i+1), x[i:(i+1)], col="red")
  }
}

#add a second loop: many runs of plotting! multiple lineages
bm.plot <- function( sigma=1, ngens=100, nlineages=100, yylim=c(-50, 50) , xxlim=c(0,100))
{
  plot(1:length(devs), devs, type = "n", col="red", 
       ylim = yylim, 
       xlim=xxlim,
       xlab="Time", ylab="Value", main="BM Simulation")
  
  for (i in 1:nlineages)       # number of lineages to simulate
  {
    x <- c(0:100)
    devs =rnorm(ngens)    # draws from a normal distribution each gen
    for (i in 1:ngens) 
    {
      x[i+1] <- x[i] + sigma*devs[i]   # BM equation   
      # step through time, increasing x a little bit each time
      lines(i:(i+1), x[i:(i+1)], col="red")   # plot line segment
    }
  }    
}

#chunk 11: now collecting the last value of each lineage
bm.plot <- function( sigma=1, ngens=100, nlineages=100, yylim=c(-50, 50) ) 
{
  plot(1:length(devs), devs, type = "n", col="red", ylim = yylim,
       xlab="Time", ylab="Value", main="BM Simulation")
  xfinal <- c(1:nlineages)  ### initialze final values
  
  for (j in 1:nlineages) # number of lineages to simulate
  {
    x <- c(0:ngens)
    devs =rnorm(ngens) # draws from a normal distribution, one per generation
    for (i in 1:ngens) 
    {
      x[i+1] <- x[i] + sigma*devs[i]     # BM equation
      # step through time, increasing x a little bit each time
      lines(i:(i+1), x[i:(i+1)], col="red")  # plot line segment
    }
    xfinal[j] <- x[ngens+1]     ### collect the last value in each lineage
  }
  return(xfinal)     ### return the vector of final values
} 

out <- bm.plot(nlineages=30)
hist(out) #probably need more values to be more normal

#next: add a lineage split. Species 1 and species 2

#source a function!
setwd("~/Projects/UW_EQG_2017/scripts/Ex3-4")
source("OU.sim.branch.R")
ls()
#view function definition
OU.sim.branch
#execute the function
OU.sim.branch()
OU.sim.branch(sigma=5)
#sigma causes a divergence in phenotype over time -- variance grows, but mean and expected value are stable

#compact version without a loop
sigma = 1
plot(1:length(devs), devs, type = "n", col="red", ylim = yylim,
     xlab="Time", ylab="Value", main="BM Simulation")
y <- c(0, cumsum(rnorm(nsteps, sd=sigma)))
lines(0:nsteps,y)

#or can use lapply instead of the for loop. lapply operates across a list
