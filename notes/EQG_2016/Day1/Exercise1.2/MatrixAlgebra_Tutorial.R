   
#Matrix algebra session

#Let’s start with the addition of two vectors.  Type or copy the text that 
#follows into the R console a line at a time, hit Enter, and examine the 
#response.
#The first and third statements define the two vectors
#Notice that R treats vectors as 1 column or 1 row matrices
a=matrix(data=c(2,3),nrow=2,ncol=1)
a
e=matrix(data=c(3,4),nrow=2,ncol=1)
e
z=a+e
z

#Now let’s try matrix addition.  We begin by defining the two matrices that we will add.
G=matrix(data=c(2,3,3,8),nrow=2,ncol=2)
G
E=matrix(data=c(1,2,2,5),nrow=2,ncol=2)
E
P=G+E
P

#Next, let’s multiply a matrix by a vector.  The following exercise computes the 
#response to directional selection on two traits, given the G-matrix and a 
#vector of selection gradients.
beta=matrix(data=c(0.2,0.4),nrow=2,ncol=1)
beta
deltazbar=G %*% beta
deltazbar

#Now, the multiplication of two matrices.  Here we depart from the sequence in 
#Arnold’s 1994 Appendix.  We will compute the inverse of a matrix and see if the 
#product of a matrix and its inverse yields the identity matrix.
Ginverse=solve(G)
Ginverse
G %*% Ginverse

#Finally, let’s solve for the eigenvalues and eigenvectors of the G-matrix.  
#What does the output mean?

eigen(G)

#To find out what the output means, let’s write some script that will 
#take a sample using a bivariate G matrix, then plot the 95% confidence 
#ellipse for that sample, and the eigenvectors of our G matrix

#First, we define our matrix and take a look at it
G=matrix(data=c(1,0.8,0.8,1),nrow=2, ncol=2,byrow=TRUE)
G

#Second, we load two libraries that we’ll need to call
library(MASS)
library(car)

#Third, we take a sample of 100 data points from the parametric version of our G matrix
data <- mvrnorm(1000, mu = c(0,0), Sigma=G)

#Fourth, we plot our data
x.range = c(-3.5,3.5)
y.range = c(-3.5,3.5) 

plot(data, xlim=x.range, ylim=y.range)

#Fifth, we plot the 95% confidence ellipse in green
null.bar=c(0,0)
ellipse(center=null.bar, shape=G, radius=2.5, center.cex=1, lwd=5, col="green", add=TRUE)

#Sixth, let’s take a look at the eigenvectors of our G matrix, so we’ll understand the script that follows
eigen(G)$vectors

#Note that the center of the plot is null.bar, the bivariate mean
null.bar


#Seventh,let’s plot the first eigenvector in blue
#We want to use the “lines” function, so we need to specify three points to connect)
n=5
delx1 = null.bar[1]  + n*eigen(G)$vectors[1,1]
delx2 = null.bar[2]  + n*eigen(G)$vectors[2,1]
delx11 = null.bar[1]  - n*eigen(G)$vectors[1,1]
delx22 = null.bar[2]  - n*eigen(G)$vectors[2,1]
x.values=c(null.bar[1] , delx1, delx11)
y.values=c(null.bar[2] , delx2, delx22)
lines(x.values, y.values, lwd=3, col="blue")

#Finally, let’s  plot the second eigenvector in red
m=5
delx1 = null.bar[1] + m*eigen(G)$vectors[1,2]
delx2 = null.bar[2] + m*eigen(G)$vectors[2,2]
delx11 =  null.bar[1] - m*eigen(G)$vectors[1,2]
delx22 =  null.bar[2]  - m*eigen(G)$vectors[2,2]
x.values=c(null.bar[1], delx1, delx11)
y.values=c(null.bar[2] , delx2, delx22)
lines(x.values, y.values, lwd=3, col="red")

#Now try changing the G matrix and see how the eigenvectors are affected

##Plot a 3 dimensional G-matrix
library(rgl)
G=matrix(data=c(1,0.9,0.8,0.9,1,0.8,0.8,0.8,1),nrow=3, ncol=3,byrow=TRUE)
G
data <- mvrnorm(1000, mu = c(0,0,0), Sigma=G)
plot3d(data)