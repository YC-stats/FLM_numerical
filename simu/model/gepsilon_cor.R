######generate stationary time series epsilon#######
gepsilon_cor<-function(n,d1,uu){##coefficient d1=0.2 or 0.5
  eta=sqrt(3)/2*rt(n,8)####rnorm(n) 4t(8)
  eps=rep(0,n)
  s=rep(0,n)
  s[1]=eta[1]
  eps[1]=0.5*s[1]*uu[1]
  for (i in 2:n){
    s[i]=d1*s[i-1]+eta[i]
    eps[i]=0.5*s[i]*uu[i]
  }
  return(eps)
}