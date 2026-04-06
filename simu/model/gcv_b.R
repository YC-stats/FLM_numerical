gcv_b <-function(y,x,P,v){
  n=length(y) #sample size
  if (n==400){
    l=seq(10^(-5),2*10^(-4),by=10^(-5))
    #change to l=seq(10^(-9),2*10^(-8),by=10^(-9)) 
    #when using legendre polynomials in MA case
  }else if (n==800){
    l=seq(10^(-9),10^(-8),by=10^(-9))
    #change to l=seq(10^(-9),10^(-8),by=10^(-9)) 
    #when using legendre polynomials in MA case
  }
  nn=length(l)
  err_squ=rep(0,nn)
  mean_squ=rep(0,nn)
  c = dim(x)[2]
  for (k in 1:nn){
    lamb=l[k]
    R=lamb*P
    M=diag(v)
    H_hat=x%*%solve(t(x)%*%M%*%x/n+R)%*%t(x)%*%M/n
    err_squ[k]=mean(((diag(1,n)-H_hat)%*%y)^2)/
      (1-mean(diag(H_hat)))^2
  }
  list(cv=err_squ,lamb=
         l[which(err_squ==min(err_squ))])
}

