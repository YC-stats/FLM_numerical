gcv_lamb_p<-function(y,x,P){
  n=length(y) #sample size
  if (n==400){
    l=seq(10^(-8),2*10^(-7),by=10^(-8))
  }else if (n==800){
    l=seq(10^(-8),10^(-7),by=10^(-8))
    #change to 10^(-9),10^(-8),by=10^(-9) when using legendre 
  }
  nn=length(l)
  err_squ=rep(0,nn)
  mean_squ=rep(0,nn)
  for (k in 1:nn){
    lamb=l[k]
    R=lamb*P
    H_hat=x%*%solve(t(x)%*%x/n+R)%*%t(x)/n
    err_squ[k]=mean(((diag(rep(1,n))-H_hat)%*%y)^2)/
      (1-mean(diag(H_hat)))^2
  }
  list(gcv=err_squ,lamb=
         l[which(err_squ==min(err_squ))])
}