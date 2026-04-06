gcv_lamb<-function(y,x,c,P){
  n=length(y) #sample size
  if (n==400){
    l=seq(10^(-8),3*10^(-7),by=10^(-8))
    #change to l=seq(10^(-9),2*10^(-8),by=10^(-9)) 
    #when using legendre polynomials in MA case
  }else if (n==800){
    l=seq(10^(-7),2*10^(-7),by=10^(-8))
    #change to l=seq(10^(-9),10^(-8),by=10^(-9)) 
    #when using legendre polynomials in MA case
  }
  nn=length(l)
  err_squ=rep(0,nn)
  #R=matrix(0,c,c)
  for (k in 1:nn){
    lamb=l[k]
    R=lamb*P
    H_hat=x[,1:c]%*%solve(t(x[,1:c])%*%x[,1:c]/n+
                R)%*%t(x[,1:c])/n
    err_squ[k]=mean(((diag(rep(1,n))-H_hat)%*%y)^2)/
                      (1-mean(diag(H_hat)))^2
    }
  list(cv=err_squ,lamb=
         l[which(err_squ==min(err_squ))])
}

