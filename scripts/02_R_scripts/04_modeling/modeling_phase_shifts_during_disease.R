for (intervals in c(seq(0,3,by=0.5)[-1])) {
 
# intervals <- 0.5 # every 2h

writeLines(paste0("resolution of data collection: ",4*intervals," h"))
  
  x<-seq(from = -3,to=3,by = intervals)
  
  plot(x,sin(x),
       main = paste0("Daily rhythms;\n", "samples collected every ", 4*intervals, " hr"),
       type = "l",
       xlab="time",
       ylab="expression",
       col="grey90", lwd=20,
       xaxt="n",
       xlim = c(-3,3),
       ylim = c(-2,2))
  axis(1, at=seq(-3,3,by=1), labels=seq(0,24,by=4))
  # axis(1, xaxp=c(-3,-1.5,0,1.5,3))
  
  
  
  for (j in c(2,3)){
    ## scenario 1: the rhythms are not inverted; show possible daily rhythms with phase shift (i)
    ## highlighted: not inverted, no phase shift
    for (i in seq(-3,3, by=1)){
      if (i %in% c(0)) {
        points(x,2*sin(x+i),type = "l",xlab="x",ylab="y",col="red", lwd=2)
      } else
        points(x,2*sin(x+i),type = "l",xlab="x",ylab="y",col="pink")  
    }
    ## scenario 2: the waveform for daily rhythms is inverted; rhythms with varying possible phase shits
    ## highlighted: inverted, and a 12h phase shift
    for(i in seq(-3,3, by=1)){
      if (i %in% c(-3,3)) {
        points(x,-0.5*sin(x+i),type = "l",xlab="x",ylab="y",col="blue", lwd=2)
      } else
        points(x,-0.5*sin(x+i),type = "l",xlab="x",ylab="y",col="lightblue")
    }
  }
}

