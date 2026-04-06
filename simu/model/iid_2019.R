source("gepsilon.R")
source("gAR0.R")
source("select_m.R")
set.seed(1687)
n=400
repli=1000
alpha=c(0.05,0.1)
res=100
a=11
d1=0.2
d2=0.5
tau=0.1
count1=0
count2=0
repli_simu=1000
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
b[1]=0.8
b[2]=0.5
b[3]=-0.3
for (i in 4:a)
  b[i]=exp(-i)
basis=matrix(1,a,res)
basis[2:a,]=sqrt(2)*cos(pi*(1:(a-1))%*%t(t))
beta_true=as.vector(t(b)%*%basis)#res dimension
err=rep(0,repli)
for (k in 1:repli){
  u=matrix(runif(n*a,-sqrt(3),sqrt(3)),n,a)
  xx=sweep(u,2,(1:a)^(-1),"*")
  #decay=exp(-0.5*(0:(a-1)))
  #xx=gAR0(n,d2,a,decay)  
  x_mat=xx%*%basis #n*res
  x_med=prcomp(x_mat,rank=a,center=TRUE)
  fpc_score=as.matrix(x_med$x/sqrt(res)) #n*rank
  fpc_basis=as.matrix(x_med$rotation*sqrt(res)) #res*rank
  fpc_eigen=apply(fpc_score^2,2,mean)
  eps=rnorm(n,0,1)
    #gepsilon(n,d1)
  y=as.vector(xx%*%b+eps)
  
  m=select_m(n,y,fpc_eigen,fpc_score)
  b_hat=apply(sweep(fpc_score[,1:m],1,y,"*"),2,mean)/fpc_eigen[1:m] #rank*1
  beta_hat=fpc_basis[,1:m]%*%b_hat #res*1
  
  sigma_hat=sqrt(mean((y-mean(y)-fpc_score[,1:m]%*%b_hat)^2))
  chi_2=matrix(rchisq(m*repli_simu,1),m,repli_simu)
  r=sqrt(apply(sweep(chi_2,1,fpc_eigen[1:m],"/"),2,sum))
  crit=sort(r)[floor(repli_simu*(c(1,1)-alpha))]
  for (i in 1:kk){
    upper[i,]=beta_hat+sigma_hat*crit[i]*sqrt(1/tau)/sqrt(n)
    lower[i,]=beta_hat-sigma_hat*crit[i]*sqrt(1/tau)/sqrt(n)
  }
  if (sum(upper[1,]>beta_true)>=(1-tau)*res&#(1-tau)*
      (sum(lower[1,]<beta_true)>=(1-tau)*res)){
    count1=count1+1
    width_valid1[k]=mean(upper[1,]-lower[1,])
  }
  if (sum(upper[2,]>beta_true)>=(1-tau)*res&
      (sum(lower[2,]<beta_true)>=(1-tau)*res)){
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

