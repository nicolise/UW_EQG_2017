times<-c(1:100)
start<-0
states<-c(start, rep(NA,length(times)))
for (time in times) {
  states[time+1]<-states[time]+RFun(1)
}

plot(c(0,times),states, type="l")
