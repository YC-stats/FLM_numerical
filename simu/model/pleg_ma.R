rm(list=ls())
source("gAR.R")
source("gMA.R")
source("gepsilon.R")
source("gepsilon_cor.R")
source("pboots1.R")
source("mv.R")
source("gcv_lamb.R")
source("leg_constru.R")
library(stats)
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
b=rep(0,a)
g=1
width_1=rep(0,repli)
width_2=rep(0,repli)
lamb=rep(0,repli)
upper=matrix(0,nrow=kk,ncol=res)
lower=matrix(0,nrow=kk,ncol=res)
b[1]=0.8
b[2]=0.5
b[3]=-0.3
for (i in 4:a)
  b[i]=i^(-3)
l=0.15 #lower multiplier for choosing window size m
#n=400 d2=0.5 (l=0.15,u=0.85) d2=1 (g=1:l=0.25,u=1.15;g=0:l=0.25,u=1.25)
#n=800 d2=0.5 (g=1:l=0.15,u=1;g=0:l=0.15,u=1.15) d2=1 (g=0:l=0.15,u=1;g=1:l=0.25,u=1.25)
u=0.85 #upper multiplier 
leg_true=matrix(0,nrow=a,ncol=res)
tilde_t=2*t-1
for (i in 1:a)
  leg_true[i,]=sqrt(2*i-1)*leg_constru(tilde_t,i-1)
beta_true=as.vector(t(leg_true)%*%b)#res dimension
for (k in 1:repli){
  cvp=0
  mm=0
  xx=gMA(n,d2,a,decay=c((1:a)^(-1)))
  x_med=xx%*%leg_true
  med=prcomp(x_med,rank=10,center=TRUE)
  u1=med$x/sqrt(res)
  uu=u1[,1]
  eps=gepsilon_cor(n,d1,uu)
  while (cvp <= 0.85){
    mm=mm+1
    cvp=summary(med)$importance[3*mm]
  }
  if (mm<=3){
    c=2*mm
  }else{
    c=2*3
  }
  y=as.vector(xx%*%b+eps)
  ####regularized version####
  f_hat=apply(xx,2,sd) #a dimension
  x=sweep(xx,2,f_hat,"/")#n*a matrix,203280
  if (c==7){
    P=diag(c(0,0,720,8400,49680,203280,655200)/(f_hat[1:c]^2),c,c)
    P[3,5]=1440*sqrt(5)/f_hat[3]/f_hat[5]
    P[6,4]=3360*sqrt(77)/f_hat[4]/f_hat[6]
    P[3,7]=1008*sqrt(65)/f_hat[3]/f_hat[7]
    P[5,7]=40320*sqrt(13)/f_hat[5]/f_hat[7]
    P[4,6]=P[6,4]
    P[5,3]=P[3,5]
    P[7,3]=P[3,7]
    P[7,5]=P[5,7]
  }else if (c==6){
    P=diag(c(0,0,720,8400,49680,203280)/(f_hat[1:c]^2),c,c)
    P[3,5]=1440*sqrt(5)/f_hat[3]/f_hat[5]
    P[6,4]=3360*sqrt(77)/f_hat[4]/f_hat[6]
    P[4,6]=P[6,4]
    P[5,3]=P[3,5]
  }else if (c==5){
    P=diag(c(0,0,720,8400,49680)/(f_hat[1:c]^2),c,c)
    P[3,5]=1440*sqrt(5)/f_hat[3]/f_hat[5]
    P[5,3]=P[3,5]
  }else if (c==4){
    P=diag(c(0,0,720,8400)/(f_hat[1:c]^2),c,c)
  }
  lamb[k]=gcv_lamb(y,x,c,P)$lamb #1.2
  theta_hat=as.vector(solve(t(x[,1:c])%*%
                              x[,1:c]/n+lamb[k]*P)%*%t(x[,1:c])%*%y/n)
  b_hat=theta_hat/f_hat[1:c]
  beta_hat=as.vector(t(leg_true[1:c,])%*%b_hat)
  ###estimated residuals
  eps_hat=as.vector(y-x[,1:c]%*%theta_hat)
  z=matrix(0,nrow=c,ncol=n)
  for (i in 1:n)
    z[,i]=t(x[i,1:c]*eps_hat[i])
  S=t(apply(z,1,cumsum))
  ####choose window size m#####
  m=mv(S,range_m=(floor(dim(S)[2]^(1/3)*l)):
         (floor(u*dim(S)[2]^(1/3))))$m
  #upper bound for m may change to 0.75(or 0.85)*original when g=0 or stronger dependence
  ####construct SCB#####
  med=pboots1(S,x=x[,1:c],f_hat=f_hat[1:c],
              t,m,c,time,lamb[k]*P,leg_true,g)
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