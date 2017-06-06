# Tutorial on direct and indirect selection
#Open R and begin the session  by specifying a P-matrix
P=matrix(data=c(1,0,0,1), nrow=2, ncol=2)
P
#We can make the specification a little less tedious by writing a function
M=function(a, b, c) {matrix(data=c(a,b,b,c),nrow=2,ncol=2)}
#so that when we call it, we do a little less typing
P=M(1,0,1)
P
#Now, let’s specify a beta vector corresponding to direct selection on both traits of equal magnitude
beta=c(1,1)
beta
#If that selection acts on our P-matrix, the resulting shifts in the means are given by the s-vector
#For this step we are using equation (2.03) from our lecture notes
s=P%*%beta
s
#If the two traits are highly correlated, we solve for the s-vector with the following steps
P=M(1,0.9,1)
P
s=P%*%beta
s
#But now if selection acts on just one trait, we get a bit of a surprise when we solve for s
beta=c(1,0)
beta
s=P%*%beta
s
#Can you draw graphs that illustrate each of these three selection scenarios?

