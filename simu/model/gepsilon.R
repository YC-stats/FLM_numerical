######generate stationary time series epsilon#######
gepsilon<-function(n,d1){##coefficient d1=0.2 or 0.5
  eta=rnorm(n,0,1)#2 for penalized version
  eps=rep(0,n)
  eps[1]=eta[1]
  for (i in 2:n)
    eps[i]=d1*eps[i-1]+eta[i]
  return(eps)#[51:n])
}