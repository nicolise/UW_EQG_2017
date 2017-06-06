###Heritability of snake vertebral numbers
###? 2011 Josef C Uyeda & Stevan J Arnold
###(revised Jan 2012)
##Stage 1_Get data into R and estimate heritability by regression

#Now, set your working directory
#by specifying the local path
setwd("~/Projects/UW_EQG_2017")
#Now read in your data
thamnophis = read.table('data/Ex1-1/R_inland_snake_data.txt',header=TRUE)

#Attach the data frame to the local environment
attach(thamnophis)

##Calculate the heritability of "body"
#Create a vector of all the unique family codes
family.list=unique(family) #161 total
#family means: from the same mother
family.list
#We will now calculate the heritability of the body vertebrae. We will need to 
#do two things, create a vector of all of the mom's vertebrae counts, and create
#a vector of all the kid's vertebrae counts averaged over families. 
#first, make sure your vector is empty:
moms=NULL
kids=NULL
#Go through a loop with i equal to a each successive value in family.list; copy from 'for' to '}' and paste in R console
for (i in family.list){ 
  byfam=subset(body,family==i) #Create a vector of body vertebrae counts for family i
  moms=rbind(moms,byfam[1]) #Save the moms count as the next element in the vector "moms"
  kids=rbind(kids,mean(byfam[-1], na.rm=TRUE)) #Save the mean of the kids count as the next element in the vector "kids" (you are excluding the 1st element of the vector, which is the mom)
}
#Linear regression of kids on moms, will save a bunch of stuff in ls; to see the summary statistics run "summary(ls)", we especially want the slope, which is stored in ls$coefficients
ls=lsfit(moms,kids) 

plot(moms, kids, col='mediumslateblue', pch=19, cex=1.5)

#Extract the 2nd coefficient from the object "coefficients" in "ls"
slope=ls$coefficients[2] 
#the heritability estimate is twice the slope because we are using just one parent, not both
h.squared=2*slope  
h.squared

#adds onto the existing plot
abline(ls, lwd=3)

##Calculate the heritability of "tail"
#Create a vector of all the unique family codes
family.list=unique(family)
#We will now calculate the heritability of the tail vertebrae. We will need to 
#do two things, create a vector of all of the mom's vertebrae counts, and create
#a vector of all the kid's vertebrae counts averaged over families. 
moms=NULL
kids=NULL
#Go through a loop with i equal to a each successive value in family.list 
for (i in family.list){ 
  byfam=subset(tail,family==i) #Create a vector of tail vertebrae counts for family i
  moms=rbind(moms,byfam[1]) #Save the moms count as the next element in the vector "moms"
  kids=rbind(kids,mean(byfam[-1], na.rm=T)) #Save the mean of the kids count as the next element in the vector "kids" (you are excluding the 1st element of the vector, which is the mom)
}
#Linear regression of kids on moms, will save a bunch of stuff in ls; to see the summary statistics run "summary(ls)", we especially want the slope, which is stored in ls$coefficients
ls=lsfit(moms, kids) 
#Extract the 2nd coefficient from the object "coefficients" in "ls"
slope=ls$coefficients[2] 
h.squared=2*slope
plot(moms, kids, col='mediumslateblue', pch=19, cex=1.5)
abline(ls, lwd=3)

#Stage 2_Use bootstrapping to make an interval estimate of heritability
##Make a function that can quickly calculate the heritability for each trait
# In the user-defined function h2, V is the name of the trait whose heritability will be estimated
h2=function(V){
  family.list=unique(family)
  moms=NULL
  kids=NULL
  for (i in family.list){
    byfam=subset(V, family==i) 
    moms=rbind(moms, byfam[1]) 
    #mean of byfam EXCEPT MOM
    kids=rbind(kids, mean(byfam[-1], na.rm=TRUE)) 
  }
  plot(moms, kids, col='mediumslateblue', pch=19, cex=1.5)
  ls=lsfit(moms, kids)
  abline(ls, lwd=3)
  
  slope=ls$coeff[2]
  h.squared=2*slope
  h.squared
}

#test new function called h2:
h2(body)
h2(ilab)
h2(tail)

##Bootstrap the heritability of "body"
#Hit 'Enter' to escape from the cycle of plots
family.list=unique(family)
h2.all=NULL
#make a new list: Samples from family.list WITH REPLACEMENT
#simulates a larger sample! can set to <1000 to run faster
#currently resamples entire families. Could also bootstrap within each level (individual daughters) then bootstrap at the higher levels (individual families). Both valid. For heritability estimate, resampling whole families is best.
for (j in 1:1000){
  moms=NULL
  kids=NULL
  family.list.boot=sample(family.list,replace=TRUE)
  for (i in family.list.boot){
    byfam=subset(body,family==i) 
    moms=rbind(moms,byfam[1]) 
    kids=rbind(kids,mean(byfam[-1],na.rm=TRUE))
  }
  
  ls=lsfit(moms, kids)
  slope=ls$coeff[2]
  h.squared=2*slope
  #Save the heritability of replicate i as the ith element in vector h2.all
  h2.all=rbind(h2.all,h.squared)  
  #Graph showing the variation in points sampled
  plot(moms,kids,xlim=c(150,180),ylim=c(150,180)) 
  #Superimpose on that graph the regression ls
  abline(ls$coefficients)
  Sys.sleep(0.1)
}
#Plot a histogram of the heritability estimates and corresponding summary stats 
hist(h2.all) 
mean(h2.all)
var(h2.all) #Our estimate of the sampling variance of the mean
sd(h2.all)  #Our estimate of the standard error of the mean

##Create a program h2boot that calculates bootstrap values so we can get the value for #each variable quickly; copy the following statements all the way thru the 3rd "}" and #paste into console
family.list=unique(family)
#R is number of boot.reps (my addition!)
h2boot=function(V,R){
  h2.all=NULL
  for (j in 1:R){
    moms=NULL
    kids=NULL
    family.list.boot=sample(family.list,replace=TRUE) #This is the key to bootstrapping
    for (i in family.list.boot){
      byfam=subset(V,family==i) 
      moms=rbind(moms,byfam[1]) 
      kids=rbind(kids,mean(byfam[-1],na.rm=TRUE))
    }
    
    ls=lsfit(moms, kids)
    slope=ls$coeff[2]
    h2=2*slope
    #Save the heritability of replicate i as the ith element in vector h2.all
    h2.all=rbind(h2.all,h2)  
    
    plot(moms,kids,xlim=c(min(V,na.rm=TRUE),max(V,na.rm=TRUE)),ylim=c(min(V,na.rm=TRUE),max(V,na.rm=TRUE))) #Graph showing the variation in points sampled
    #Superimpose on that graph the regression ls
    abline(ls$coefficients)
    Sys.sleep(0.1) #tells the systemto pause on each plot to allow user to visualize
  }
  hist(h2.all) 
  #Use a double-headed arrow so that out is a global variable not just a local variable
  out<<-list(mean=mean(h2.all),sd=sd(h2.all)) 
  out
}

#Use 'h2boot(tail,1000)' to start the program on 'tail' with 1000 boot.reps; hit 'ESC' to end the cycle of boot plots
#but if you hit 'ESC' the histogram and stats won't appear
h2boot(tail,1000)
h2boot(ilab,100)
