######generate stationary FMA(1)#######
gMA<-function(n,d2,a,decay){
  xx=matrix(0,nrow=a,ncol=n) #nrow=(a+50)
  ee=matrix(rnorm(a*n),nrow=a,ncol=n)
  e=sweep(ee,1,decay,'*') #column
  AA=diag(a)
  AA[abs(row(AA)-col(AA))==1]=1/5
  eps=AA%*%e
  xx[,1]=eps[,1]
  for (i in 2:n) #(a+50)
    xx[,i]=eps[,i]+d2%*%e[,(i-1)]
  return(t(xx))
}