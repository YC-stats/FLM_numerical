rm(list=ls())
source("gMA.R")
source("gepsilon0.R")
source("pboots_p.R")
source("mv.R")
source("gcv_lamb_p.R")
#source("leg_constru.R")
library(matrixcalc)
set.seed(1687)
n=400
repli=1000
time=1000
d2=0.2 #dependence of AR or MA process
alpha=c(0.05,0.1)
res=100
count_1=0
count_2=0
p=1
r=matrix(0,time,p)
width_1=matrix(0,repli,p)
width_2=matrix(0,repli,p)
t=(1:res)/res
kk=length(alpha)
a=10
c=rep(0,p)
b=rep(0,a)
upper=array(0,c(kk,res,p))
lower=array(0,c(kk,res,p))
rank1=6
g_hat=matrix(0,res,p)
xx=array(0,c(n,a,p))
lamb=rep(0,repli)
basis=matrix(1,a,res)
basis[2:a,]=sqrt(2)*cos(pi*(1:(a-1))%*%t(t))
g=0
beta_true=cbind()
b_comb=cbind()
beta_hat=matrix(0,res,p)
b[1]=0.8
b[2]=0.5
b[3]=-0.3
for (ss in 4:a)
  b[ss]=ss^(-4)
beta=as.vector(t(b)%*%basis) #res*1 vector
for (i in 1:p){
  beta_true=cbind(beta_true,beta) #res*p
  b_comb=c(b_comb,b) #a*p
}
l=0.15 #lower multiplier for choosing window size m
#n=800 p=3 (l=0.25,u=1.15) p=2 (l=0.65,u=1.15) p=1 (l=0.25,u=1)
#n=400 p=3 (l=0.25,u=1.15) p=2 (l=0.15,u=0.75) p=1 (l=0.25,u=1.15)
u=1.15 #upper multiplier 
for (k in 1:repli){
  f=cbind()
  P_comb=list()
  xx_comb=cbind()
  x_comb=cbind() #scaled version
  for (i in 1:p){
    cvp=0
    mm=0
    xx=gMA(n,d2,a,decay=(1:a)^(-1.2)) #n*a
    x_mat=xx%*%basis #n*res
    x_med=prcomp(x_mat,rank=rank1,center=TRUE)
    while (cvp <= 0.85){
      mm=mm+1
      cvp=summary(x_med)$importance[3*mm]
    }
    c[i]=2*mm
    f_hat=as.vector(apply(xx[,1:c[i]],2,sd))
    pen_mat=list(diag((c(0,2*(pi*(1:(c[i]-1)))^4))/
                           (f_hat^2),c[i],c[i]))
    P_comb=append(P_comb,pen_mat)
    x_med=sweep(xx[,1:c[i]],2,f_hat,"/")
    xx_comb=cbind(xx_comb,xx)
    x_comb=cbind(x_comb,x_med)# scaled version
    f=c(f,f_hat)#vector pc*1 #combined
  }
  eps=rnorm(n,0,1)#gepsilon0(n,d1)
  y=as.vector(xx_comb%*%b_comb+eps)
  x=x_comb #n*(c[1]+...+c[p])
  c_length=sum(c)
  pen_matrix=matrix(0,c_length,c_length)
  if (p==1){
    pen_matrix[1:c[1],1:c[1]]=P_comb[[1]]
  }else if (p==2){
    pen_matrix[1:c[1],1:c[1]]=P_comb[[1]]
    pen_matrix[(c[1]+1):(c[1]+c[2]),(c[1]+1):(c[1]+c[2])]=P_comb[[2]]
  }else{
    pen_matrix[1:c[1],1:c[1]]=P_comb[[1]]
    pen_matrix[(c[1]+1):(c[1]+c[2]),(c[1]+1):(c[1]+c[2])]=P_comb[[2]]
    pen_matrix[(c[1]+c[2]+1):(c[1]+c[2]+c[3]),
               (c[1]+c[2]+1):(c[1]+c[2]+c[3])]=P_comb[[3]]
  }
  lamb[k]=gcv_lamb_p(y,x,pen_matrix)$lamb
  P=lamb[k]*pen_matrix
  theta_hat=as.vector(solve(t(x)%*%x/n+P)%*%t(x)%*%y/n)
  b_hat=as.vector(theta_hat/f)
  if (p==1){
    beta_hat[,1]=as.vector(t(b_hat)%*%basis[1:c[1],])
  }else if (p==2){
    beta_hat[,1]=as.vector(t(b_hat[1:c[1]])%*%basis[1:c[1],])
    beta_hat[,2]=as.vector(t(b_hat[(c[1]+1):c_length])%*%basis[1:c[2],])
  }else{
    beta_hat[,1]=as.vector(t(b_hat[1:c[1]])%*%basis[1:c[1],])
    beta_hat[,2]=as.vector(t(b_hat[(c[1]+1):(c[1]+c[2])])%*%basis[1:c[2],])
    beta_hat[,3]=as.vector(t(b_hat[(c[1]+c[2]+1):c_length])%*%basis[1:c[3],])
    }
  eps=as.vector(y-x%*%theta_hat)
  z=matrix(0,c_length,n)
  for (i in 1:n)
    z[,i]=as.vector(x[i,]*eps[i])
  S=t(apply(z,1,cumsum))##(5a)*n matrix
  ####choose window size
  m=mv(S,range_m=(floor(dim(S)[2]^(1/3)*l)):
         (floor(u*dim(S)[2]^(1/3))))$m
  ####construct SCB#####
  med=pboots_p(S,x,f,t,m,c,time,P,basis,g,p)
  count1=0
  count2=0
  if (g==1){
    r=sort(med)
    crit=r[floor(time*(c(1,1)-alpha))]
    for (i in 1:p){
      for (j in 1:kk){
        upper[j,,i]=beta_hat[,i]+crit[j]/sqrt(n)#crit[i]
        lower[j,,i]=beta_hat[,i]-crit[j]/sqrt(n)
      }
      if (sum(upper[1,,i]>beta_true[,i])==res&
          (sum(lower[1,,i]<beta_true[,i])==res)){
        count1=count1+1
        width_1[k]=mean(upper[1,,i]-lower[1,,i])
      }
      if (sum(upper[2,,i]>beta_true[,i])==res&
          (sum(lower[2,,i]<beta_true[,i])==res)){
        count2=count2+1
        width_2[k]=mean(upper[2,,i]-lower[2,,i])
      }
    }
    if (count1==p) count_1=count_1+1
    if (count2==p) count_2=count_2+1
  }else{
    r=sort(med$M_r)
    crit=r[floor(time*(c(1,1)-alpha))]
    g_hat=med$g_hat
    for (i in 1:p){
      for (j in 1:kk){
        upper[j,,i]=beta_hat[,i]+crit[j]*g_hat[,i]/sqrt(n)#crit[i]
        lower[j,,i]=beta_hat[,i]-crit[j]*g_hat[,i]/sqrt(n)
      }
      if (sum(upper[1,,i]>beta_true[,i])==res&
          (sum(lower[1,,i]<beta_true[,i])==res)){
        count1=count1+1
      }
      width_1[k,i]=mean(upper[1,,i]-lower[1,,i])
      if (sum(upper[2,,i]>beta_true[,i])==res&
          (sum(lower[2,,i]<beta_true[,i])==res)){
        count2=count2+1
      }
      width_2[k,i]=mean(upper[2,,i]-lower[2,,i])
    }
    if (count1==p) count_1=count_1+1
    if (count2==p) count_2=count_2+1
  }
}
print(count_1/repli)
print(count_2/repli)
width1=width_1[,1]
width2=width_2[,1]
print(round(mean(width1[width1>0]),4))
print(round(mean(width2[width2>0]),4))
