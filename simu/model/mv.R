######minimum volatility method to choose m#######
mv<-function(S,range_m){
  r=dim(S)[1]#p*a
  n=dim(S)[2]#n
  mm=length(range_m)
  SE=rep(0,(mm-4))
  Xi=array(0,dim=c(r,r,mm))
  #SS r*(n-m+1) matrix
  for (i in 1:mm){
    m=range_m[i]
    SS=cbind(S[,m],S[,(m+1):n]-S[,1:(n-m)])/sqrt(m)
    Xi[,,i]=SS%*%t(SS)/(n-m+1)
  }
  Xi_bar=array(0,dim=c(r,r,(mm-4)))
  for (i in 3:(mm-2)){
    Xi_bar[,,(i-2)]=
      apply(Xi[,,(i-2):(i+2)],c(1,2),mean)
    sqsum=0
    for (j in (i-2):(i+2))
      sqsum=sqsum+sum((Xi[,,j]-Xi_bar[,,(i-2)])^2)
    SE[i-2]=(sqsum/4)^(1/2)
  }
  list(vseries=SE,m=range_m[which(SE==min(SE))+2])
}