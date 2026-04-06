######generate stationary time series epsilon#######
gepsilon0<-function(n,d1){
  eta=rnorm(n,0,1)
  eps=rep(0,n)
  eps[1]=eta[1]
  for (i in 2:n)
    eps[i]=d1*eps[i-1]+eta[i]
  return(eps)#[51:n])
}