### Several Simulation functions
#   bm.plot = basic bm simulation plot
#   bm.plot.fast = faster basic bm sim
#   ou.plot = basic ou sim
#   OU.sim.branch = OU simulation with branch and different optima
#                   sourced from OU.sim.branch.R
#################################################

bm.plot <- function( sigma=1, ngens=100, nlineages=100, yylim=c(-50, 50) ) 
{
    plot(1:ngens, 1:ngens, type = "n", col="red", ylim = yylim,
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
    return(xfinal)     ### return the final values
} 


bm.plot.fast <- function( sigma=1, nsteps=100, nlineages=100, yylim=c(-50, 50))
{
# Set up plotting environment
	plot(0, 0, type = "n", xlab = "Time", ylab = "Trait", 
		xlim=c(0, nsteps), ylim=yylim)

# Draw random deviates and plot
	lapply( 1:nlineages, function(x) 
	  lines(0:nsteps, c(0, cumsum(rnorm(nsteps, sd=sigma)))))
}

source("OU.sim.branch.R")
