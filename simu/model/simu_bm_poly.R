source("pboots1.R") 
source("mv.R")
source("gcv_lamb.R")
source("poly_constru.R")
source("gepsilon.R")
source("basis_poly.R")
library(MASS)
set.seed(1687)
n=400
repli=1000
time=1000
alpha=c(0.05,0.1)
res=100
a=10
d1=0.2
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
width_valid1=rep(0,repli)
width_valid2=rep(0,repli)
g=0
rank1=10
b[1]=0.8
b[2]=0.5
b[3]=-0.3
for (i in 4:a)
  b[i]=exp(-i)
basis_type="Polynomial" #change to Legen, FPC
if (basis_type=="Four"){
  basis=matrix(0,a,res)
  basis[1,]=sqrt(3)*t
  basis[2:a,]=sqrt(2)*sin(pi*(1:(a-1))%*%t(t))
}else if (basis_type=="Polynomial"){
  basis=matrix(0,a,res)
  for (i in 1:a)
    basis[i,]=poly_constru(t,i-1)
  }
beta_true=as.vector(t(b)%*%basis) #res dimension
beta_hat=matrix(0,res,repli)
err=rep(0,repli)
scaled_vec=rep(0,a)
for (k in 1:a)
  scaled_vec[k]=1/sqrt(k*(k+1))
for (k in 1:repli){
  cvp=0
  mm=0
  gbasis=matrix(0,a,res)
  for (i in 1:a)
    gbasis[i,]=basis_poly(t,i-1)
  if (basis_type=="Four"){
    xx=sweep(matrix(rnorm(n*a,0,1),n,a),2,c(1,pi*(1:(a-1))),"/")
    x_mat=t(sapply(1:n,function(i){
      xx[i,1]*t+apply(sweep(basis[2:a,],
                            1,xx[i,2:a],"*"),2,sum)},simplify = "cbind"))
  }else if (basis_type=="Polynomial"){
    x_score=sweep(matrix(rnorm(n*a,0,1),n,a),2,scaled_vec,"*")
    x_mat=sapply(1:n,function(i){
      apply(sweep(gbasis,1,x_score[i,],"*"),2,sum)},
      simplify="rbind")
  }
  x_med=prcomp(x_mat,rank=10,center=TRUE)
  while (cvp <= 0.85){
    mm=mm+1
    cvp=summary(x_med)$importance[3*mm]
  }
  c=2*mm
  eps=#rnorm(n,0,0.25)
    gepsilon(n,d1)
  xx=t(ginv(t(basis))%*%x_mat)
  y=as.vector(xx%*%b+eps)
  ####regularized version####
  f_hat=apply(xx,2,sd) 
  x=sweep(xx,2,f_hat,"/")
  if (basis_type=="Four"){
    P=diag((c(0,(pi*(1:(c-1)))^4))/(f_hat[1:c]^2),c,c)
  }else if (basis_type=="Polynomial"){
    if (c==8){
      P=diag(c(0,120,2520,15120,55440,154440,360360,742560)/(f_hat[1:c]^2),c,c)
      P[2,8]=1783.928/f_hat[2]/f_hat[8]
      P[4,8]=49049.92/f_hat[4]/f_hat[8]
      P[6,8]=247863.3/f_hat[6]/f_hat[8]
      P[2,4]=549.9091/f_hat[2]/f_hat[4]
      P[2,6]=1111.539/f_hat[2]/f_hat[6]
      P[3,5]=6473.997/f_hat[3]/f_hat[5]
      P[3,7]=11389.89/f_hat[3]/f_hat[7]
      P[4,6]=30562.28/f_hat[4]/f_hat[6]
      P[5,7]=97537.19/f_hat[5]/f_hat[7]
      P[8,2]=P[2,8]
      P[8,4]=P[4,8]
      P[8,6]=P[6,8]
      P[4,2]=P[2,4]
      P[6,2]=P[2,6]
      P[5,3]=P[3,5]
      P[7,3]=P[3,7]
      P[6,4]=P[4,6]
      P[7,5]=P[5,7]
    }else if (c==6){
      P=diag(c(0,120,2520,15120,55440,154440)/(f_hat[1:c]^2),c,c)
      P[2,4]=549.9091/f_hat[2]/f_hat[4]
      P[2,6]=1111.539/f_hat[2]/f_hat[6]
      P[3,5]=6473.997/f_hat[3]/f_hat[5]
      P[4,6]=30562.28/f_hat[4]/f_hat[6]
      P[4,2]=P[2,4]
      P[6,2]=P[2,6]
      P[5,3]=P[3,5]
      P[6,4]=P[4,6]
    }else if (c==4){
      P=diag(c(0,120,2520,15120)/(f_hat[1:c]^2),c,c)
      P[2,4]=549.9091/f_hat[2]/f_hat[4]
      P[4,2]=P[2,4]
    }
  }
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
  m=mv(S,range_m=c(1:
                     floor(0.75*dim(S)[2]^(1/3))))$m
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
    width_valid1[k]=mean(upper[1,]-lower[1,])
  }
  if (sum(upper[2,]>beta_true)==res&
      (sum(lower[2,]<beta_true)==res)){
    count2=count2+1
    width_valid2[k]=mean(upper[2,]-lower[2,])
  }
  width_1[k]=mean(upper[1,]-lower[1,])
  width_2[k]=mean(upper[2,]-lower[2,])
}
print(count1/repli)
print(count2/repli)
print(round(mean(width_1),3))
print(round(mean(width_2),3))

