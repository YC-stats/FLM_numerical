gcv_lamb<-function(y,x,c,P){
  n=length(y) #sample size
  l=seq(10^(-13),10^(-10),by=10^(-13)) #10^(-13) to 10^(-10) for Legendre  10/8 for FPC
  nn=length(l)
  err_squ=rep(0,nn)
  R=matrix(0,c,c)
  for (k in 1:nn){
    lamb=l[k]
    R=lamb*P
    H_hat=x%*%solve(t(x)%*%x/n+R)%*%t(x)/n
    err_squ[k]=mean(((diag(rep(1,n))-H_hat)%*%y)^2)/
      (1-mean(diag(H_hat)))^2
  }
  list(cv=err_squ,lamb=
         l[which(err_squ==min(err_squ))])
}

