rm(list=ls())
source("gAR.R")
source("gMA.R")
source("gepsilon_cor.R")
#source("gepsilon.R")
source("pboots1.R") 
source("mv.R")
source("gcv_lamb.R")
set.seed(1687)#change to 687 when n=800
n=400
repli=1000
time=1000
d1=0.2  
d2=0.5
alpha=c(0.05,0.1)
res=100
a=10
count1=0
count2=0
t=(1:res)/res
kk=length(alpha)
upper=matrix(0,nrow=kk,ncol=res)
lower=matrix(0,nrow=kk,ncol=res)
b=rep(0,a)
lamb=rep(0,repli)
width_1=rep(0,repli)
width_2=rep(0,repli)
rank1=10
b[1]=0.8
b[2]=0.5
b[3]=-0.3
for (i in 4:a)
  b[i]=i^(-3)#exp(-i)
basis=matrix(1,a,res)
basis[2:a,]=sqrt(2)*cos(pi*(1:(a-1))%*%t(t))
#basis=cos(2*pi*(1:a)%*%t(t))+sin(2*pi*(1:a)%*%t(t))
beta_true=as.vector(t(b)%*%basis)#res dimension
beta_hat=matrix(0,res,repli)
err=rep(0,repli)
g=1
l=0.25 #lower multiplier for choosing window size m
#n=400 d2=0.5 (g=0:l=0.15,u=0.75; g=1:l=0.25,u=1) d2=1 (g=1:l=0.25,u=1.25;g=0:l=0.25,u=1)
#n=800 d2=0.5 (g=0:l=0.75,u=1.15;g=1:l=0.75,u=1.5) d2=1 (g=0:l=0.5,u=1.25,g=1l=0.75,u=1.5)
u=1 #upper multiplier 
for (k in 1:repli){
  cvp=0
  mm=0
  xx=gMA(n,d2,a,decay=c((1:a)^(-1)))
  x_mat=xx%*%basis
  x_med=prcomp(x_mat,rank=10,center=TRUE)
  while (cvp <= 0.85){
    mm=mm+1
    cvp=summary(x_med)$importance[3*mm]
  }
  c=2*mm
  u1=as.matrix(x_med$x/sqrt(res))
  uu=u1[,1]
  eps=gepsilon_cor(n,d1,uu)
  y=as.vector(xx%*%b+eps)
  ####regularized version####
  f_hat=apply(xx,2,sd) 
  x=sweep(xx,2,f_hat,"/")
  P=diag((c(0,2*(pi*(1:(c-1)))^4))/(f_hat[1:c]^2),c,c)
  lamb[k]=gcv_lamb(y,x,c,P)$lamb
  theta_hat=as.vector(solve(t(x[,1:c])%*%
              x[,1:c]/n+lamb[k]*P)%*%t(x[,1:c])%*%y/n)
  b_hat=theta_hat/f_hat[1:c]
  ###estimated functional coefficient
  beta_hat=as.vector(t(b_hat)%*%basis[1:c,])
  ##estimated residuals##
  eps_hat=as.vector(y-x[,1:c]%*%theta_hat)
  z=matrix(0,nrow=c,ncol=n)
  for (i in 1:n)
    z[,i]=t(x[i,1:c]*eps_hat[i])
  S=t(apply(z,1,cumsum))##a_est*n matrix
  ####choose window size m#####
  m=mv(S,range_m=c(floor(dim(S)[2]^(1/3)*l)):
         (floor(u*dim(S)[2]^(1/3))))$m
  ####construct SCB#####
  med=pboots1(S,x=x[,1:c],f_hat=f_hat[1:c],
              t,m,c,time,lamb[k]*P,basis,g)
  if (g==1){
    r=sort(med)
    crit=r[floor(time*(c(1,1)-alpha))]
    for (i in 1:kk){
      upper[i,]=beta_hat+crit[i]/sqrt(n)
      lower[i,]=beta_hat-crit[i]/sqrt(n)
    }
  }else{
    r=sort(med$M_r)
    g_hat=med$g_hat
    crit=r[floor(time*(c(1,1)-alpha))]
    for (i in 1:kk){
      upper[i,]=beta_hat+g_hat*crit[i]/sqrt(n)
      lower[i,]=beta_hat-g_hat*crit[i]/sqrt(n)
    }
  }
  if (sum(upper[1,]>beta_true)==res&
      (sum(lower[1,]<beta_true)==res)){
    count1=count1+1
  }
  if (sum(upper[2,]>beta_true)==res&
      (sum(lower[2,]<beta_true)==res)){
    count2=count2+1
  }
  width_1[k]=mean(upper[1,]-lower[1,])
  width_2[k]=mean(upper[2,]-lower[2,])
}
print(count1/repli)
print(count2/repli)
print(round(mean(width_1),3))
print(round(mean(width_2),3))

