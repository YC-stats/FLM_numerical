######generate stationary FAR(1)#######
gAR0<-function(n,d2,a,decay){
  xx=matrix(0,nrow=a,ncol=n)#nrow=(a+50)
  ee=matrix(rnorm(a*n,0,1),nrow=a,ncol=n)
  e=sweep(ee,1,decay,'*')#by row
  xx[,1]=e[,1]
  AA=diag(d2,a,a)
  for (i in 2:n)#(a+50)
    xx[,i]=e[,i]+AA%*%xx[,(i-1)]
  return(t(xx))#[,51:n])) #n*res
}