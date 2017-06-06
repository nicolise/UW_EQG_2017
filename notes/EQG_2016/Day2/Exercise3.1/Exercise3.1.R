##Computational Exercise 3.1: Estimating and plotting a selection surface
##Copyright 2016 Stevan J. Arnold

##Get your data into R using the following steps.
#Save the file radix5.txt to your desktop or folder of choice
#In the following example, we have created a folder called ?R wd? and placed the data file in it
#Now, set your working directory using a statement comparable to this one
#Be sure and use the / not the \ spacing convention!
setwd('C:/Documents and Settings/arnolds/Desktop/R wd')

#Now read in your data
thamnophis = read.table('radix5.txt',header=TRUE)
#Attach the data frame to the local environment
attach(thamnophis)
#Print out the data frame to make sure it looks OK
thamnophis

##Stage 1  Estimating the selection gradients
#Use the following statements to create standardized versions of BODY, TAIL, and TIME
#We want to standardize BODY and TAIL so that their means are zero and their sds are 1
#We want to standardize SPEED so that it's mean is 1
#These standardizations will simplify the interpretations of the coefficients that we will estimate
new.body=(BODY-mean(BODY))/sd(BODY)
new.tail=(TAIL-mean(TAIL))/sd(TAIL)
new.speed=SPEED/mean(SPEED)
#Do these transformations achieve the right means and standard deviations for the new variables?
#Find out by using statements like mean(new.body), sd(new.body).

#Let us begin by fitting a plane to the transformed data
model <-lm(new.speed~new.body + new.tail)
#print out the coefficients and statistics of the fit
#this output will give us our best estimates of the directional selection gradients for 
#new.body and new.tail
summary(model)

#Now, fit a full quadratic model to the data, remembering to use the factor of 0.5 for the 
#stabilizing selection gradients
#This model will give us estimates of the nonlinear selection gradients, gamma
#Create a new variable called prod which is the product of new.body and new.tail
prod = new.body*new.tail
#Then fit the quadratic model
model <- lm(new.speed ~ new.body + new.tail + I(0.5*new.body^2) + I(0.5*new.tail^2) + prod )
summary(model)

#In summary, using z1 for body and z2 for tail,
#Using results for beta from the linear fit and the results for gamma from the quadratic fit, we have
beta1= 0.031408  
beta2= 0.001504
gamma11 = -0.010578 
gamma22 = -0.005914 
gamma12 =  0.079017
#Because of our standardization, beta1 tells us that if we increase the value of new.body by 1 sd, 
#we will increase new.speed by 3.1%
   
##Stage 2  Plotting the selection surface, an ISS
#Set the parameters for a full quadratic fitness surface for two traits
z1 = new.body
z2 = new.tail
#Set the value for the intercept of the surface when z1=z2=0
alpha=1
#Define a function called fit that will compute the value of fitness as a function of z1 and z2
fit <-  function(z1,z2, alpha, beta1, beta2, gamma11, gamma22, gamma12)  alpha + (beta1*z1) + (beta2*z2) + (gamma11*0.5*(z1^2)) + (gamma22*0.5*(z2^2)) + (gamma12*(z1 * z2)) 
#Define a series of values for z1 and z2 that will be used to compute the value of relative fitness on the surface
x <- seq(-2, 2, length = 30)
y <- seq(-2, 2, length = 30)

# Compute the surface values of relative fitness using the x-y grid of values for z1 and 
# z2 for later use by the surface plotting function called persp
#We could use the fit function to compute values of fitness but instead we will use a #function, called outer, that is more compatible with persp
# The function outer has some specific requirements so we oblige by writing our fitness 
# function in the following form

z <- outer(x, y, function(a, b, alpha, beta1, beta2, gamma11, gamma22, gamma12) 1 + (0.031408*a) + (0.001504*b) + (-0.010578 *0.5*(a^2)) + (-0.005914 *0.5*(b^2)) + ( 0.079017 *(a * b)))

#Define two variables that give the number of rows and cols in z
nrz <- nrow(z)
ncz <- ncol(z)

# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("yellow", "orange") )
 
# Generate the desired number of colors from this palette
nbcol <- 100
color <- jet.colors(nbcol)

# Compute the z-value at the facet centres
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]

# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)

#Finally, plot the ISS
par(bg = "white")
persp(x, y, z, col=color[facetcol], phi=30, theta=-30)

##Stage 3 ? Can you build a loop that will show the contour plot while bootstrapping 
#over the snake sample? Here?s an example of script that does that!
#We begin by writing a function that we will use during bootstrapping to estimate the coefficients that describe the surface for a particular boot sample
est.coeff=function(y,x1,x2){
  prod=x1*x2
  m1=lm(y~x1+x2)
  m2=lm(y~x1+x2+I(0.5*x1^2)+I(0.5*x2^2)+prod)
  c(m1$coeff[2],m1$coeff[3],m2$coeff[4],m2$coeff[5],m2$coeff[6])
}
  

##Next, we program a perspective plot function
x <- seq(-2, 2, length = 30)
y <- seq(-2, 2, length = 30)
# Compute the surface values of relative fitness using the x-y grid of values for z1 and 
# z2 for later use by the surface plotting function called persp
#We could use the fit function to compute values of fitness but instead we will use a #function, called outer, that is more compatible with persp
# The function outer has some specific requirements so we oblige by writing our fitness 
# function in the following form
##Write a function that plots the perspective plot for each bootstrap
persp.boot=function(b1,b2,y1,y2,y12,x=seq(-10,10,length=30),y=seq(-10,10,length=30)){
z <- outer(x, y, function(a, b, alpha, beta1, beta2, gamma11, gamma22, gamma12) 1 + ( b1*a) + (b2*b) + (y1*0.5*(a^2)) + (y2*0.5*(b^2)) + (y12*(a * b)))
nrz <- nrow(z)
ncz <- ncol(z)
jet.colors <- colorRampPalette( c("yellow", "orange") )
nbcol <- 100
color <- jet.colors(nbcol)
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)
par(bg = "white")
persp(x, y, z, xlab="Number of body vertebrae", ylab="No. tail vertebrae", zlab="Crawling speed", col=color[facetcol],phi=30,theta=-30)
}

##Test it out
persp.boot(-2,2,1,1,1)

##Bootstrap over individuals
par(ask=FALSE)
n=length(new.body)
for (i in 1:100){
  samp=sample(1:n,n,replace=TRUE)
  boot.speed=new.speed[samp]
  boot.body=new.body[samp]
  boot.tail=new.tail[samp]
  boot.coeff=est.coeff(boot.speed,boot.body,boot.tail)
  persp.boot(boot.coeff[1],boot.coeff[2],boot.coeff[3],boot.coeff[4],boot.coeff[5],x=seq(-10,10,length=10),y=seq(-10,10,length=10))
Sys.sleep(0.1)
}

## Stage 4 Plotting eigenvectors on the surface
#For another view of the surface, do a simple contour plot

par(pty ="s")
contour(x, y, z)

#Let?s plot the eigenvectors of the gamma-matrix on this contour plot
#First, let?s define the matrix
gamma = matrix(c(gamma11, gamma12, gamma12, gamma22), nrow=2, ncol=2)
gamma
eigen(gamma)
#Next, let?s setup the beta vector and check to see if its values are correct
beta = matrix(c(beta1, beta2), nrow=2, ncol=1)
beta
#Using gamma and beta, we can solve for the stationary point on the surface using
#expression(11) from Phillips & Arnold 1989
z.zero = -(solve(gamma)) %*% beta
#Plot the stationary point on the surface
points(z.zero[1,1], z.zero[2,1])
#Finally, plot the eigenvectors on the surface
#First one eigenvector
delx1 = z.zero[1] + 6*eigen(gamma)$vectors[1,2]
delx2 = z.zero[2] + 6*eigen(gamma)$vectors[2,2]
delx11 =  z.zero[1] - 6*eigen(gamma)$vectors[1,2]
delx22 =  z.zero[2]  - 6*eigen(gamma)$vectors[2,2]
x.values=c(delx1, z.zero[1], delx11)
y.values=c(delx2, z.zero[2], delx22)
lines(x.values, y.values, lty=2,lwd=2)
#then the other
delx1 = z.zero[1]  + 6*eigen(gamma)$vectors[1,1]
delx2 = z.zero[2]  + 6*eigen(gamma)$vectors[2,1]
delx11 = z.zero[1]  - 6*eigen(gamma)$vectors[1,1]
delx22 = z.zero[2]  - 6*eigen(gamma)$vectors[2,1]
x.values=c (delx1, z.zero[1], delx11)
y.values=c(delx2, z.zero[2], delx22)
lines(x.values, y.values, lty=2, lwd=2)
#Which of these eigenvectors is gamma max?
#Now, let?s estimate the omega-matrix, approximating it as the negative inverse of the 
#gamma-matrix
omega=-(solve(gamma))
omega
#It?s eigenvectors should be the same as those of the  negative gamma-matrix
#Let?s check
eigen(omega)
#What is the interpretation of the eigenvalues on the contour plot?
