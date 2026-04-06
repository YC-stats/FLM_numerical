source("select_m_comb.R")
time=10000
alpha=c(0.05,0.1)
kk=length(alpha)
repli_simu=10000
tau=0.1
xx1=as.matrix(read.csv("curve_demand_2015.csv", header=FALSE))
xx2=as.matrix(read.csv("curve_demand_2016.csv", header=FALSE))
xx3=as.matrix(read.csv("curve_demand_2017.csv", header=FALSE))
y1=as.matrix(read.csv("mean_demand_2015.csv", header=FALSE))
y2=as.matrix(read.csv("mean_demand_2016.csv", header=FALSE))
y3=as.matrix(read.csv("mean_demand_2017.csv", header=FALSE))
xx=cbind(xx1,xx2,xx3)
n=dim(xx)[2] #782
y_med=as.vector(rbind(y1,y2,y3))
xx=xx[,-n]
y_ori=log(y_med[4:n])
y=as.vector(y_ori-mean(y_ori))
nn=dim(xx)[2]
res=dim(xx)[1]
upper1=matrix(0,kk,res)
lower1=matrix(0,kk,res)
upper2=matrix(0,kk,res)
lower2=matrix(0,kk,res)
upper3=matrix(0,kk,res)
lower3=matrix(0,kk,res)
temp_1=xx[,-c(1,2)] #X_{i1}(t) #res*n1
temp_2=xx[,-c(1,nn)]#X_{i2}(t)
temp_3=xx[,-c((nn-1),nn)]#X_{i3}(t)
mm=dim(temp_1)[2]
xx_1=matrix(0,res,mm)
xx_2=matrix(0,res,mm)
xx_3=matrix(0,res,mm)
for (i in 1:res){#centering
  xx_1[i,]=temp_1[i,]-apply(temp_1,1,mean)[i]#res*n matrix
  xx_2[i,]=temp_2[i,]-apply(temp_2,1,mean)[i]
  xx_3[i,]=temp_3[i,]-apply(temp_3,1,mean)[i]
}
x_med1=prcomp(t(xx_1),rank=10,center=TRUE)
x_med2=prcomp(t(xx_2),rank=10,center=TRUE)
x_med3=prcomp(t(xx_3),rank=10,center=TRUE)
fpc_score1=as.matrix(x_med1$x/sqrt(res)) #n*rank
fpc_score2=as.matrix(x_med2$x/sqrt(res)) #n*rank
fpc_score3=as.matrix(x_med3$x/sqrt(res)) #n*rank
fpc_basis1=as.matrix(x_med1$rotation*sqrt(res)) #res*rank
fpc_basis2=as.matrix(x_med2$rotation*sqrt(res))
fpc_basis3=as.matrix(x_med3$rotation*sqrt(res))
fpc_eigen1=apply(fpc_score1^2,2,mean)
fpc_eigen2=apply(fpc_score2^2,2,mean)
fpc_eigen3=apply(fpc_score3^2,2,mean)

m=max(which.min(select_m_comb(n,y,fpc_eigen1,fpc_score1)+
                  select_m_comb(n,y,fpc_eigen2,fpc_score2)+
                  select_m_comb(n,y,fpc_eigen3,fpc_score3)),2)
b_hat1=apply(sweep(fpc_score1[,1:m],1,y,"*"),2,mean)/fpc_eigen1[1:m] #rank*1
b_hat2=apply(sweep(fpc_score2[,1:m],1,y,"*"),2,mean)/fpc_eigen2[1:m]
b_hat3=apply(sweep(fpc_score3[,1:m],1,y,"*"),2,mean)/fpc_eigen3[1:m]
beta_hat1=fpc_basis1[,1:m]%*%b_hat1 #res*1
beta_hat2=fpc_basis2[,1:m]%*%b_hat2
beta_hat3=fpc_basis3[,1:m]%*%b_hat3

sigma_hat=sqrt(mean((y-mean(y)-fpc_score1[,1:m]%*%b_hat1-
                       fpc_score2[,1:m]%*%b_hat2-
                       fpc_score3[,1:m]%*%b_hat3)^2))
chi_2=matrix(rchisq(3*m*repli_simu,1),3*m,repli_simu)
r=sqrt(apply(sweep(chi_2,1,c(fpc_eigen1[1:m],fpc_eigen2[1:m],
                             fpc_eigen3[1:m]),"/"),2,sum))
crit=sort(r)[floor(repli_simu*(c(1,1)-alpha))]
for (i in 1:kk){
  upper1[i,]=beta_hat1+sigma_hat*crit[i]*sqrt(1/tau)/sqrt(n)
  lower1[i,]=beta_hat1-sigma_hat*crit[i]*sqrt(1/tau)/sqrt(n)
  upper2[i,]=beta_hat2+sigma_hat*crit[i]*sqrt(1/tau)/sqrt(n)
  lower2[i,]=beta_hat2-sigma_hat*crit[i]*sqrt(1/tau)/sqrt(n)
  upper3[i,]=beta_hat3+sigma_hat*crit[i]*sqrt(1/tau)/sqrt(n)
  lower3[i,]=beta_hat3-sigma_hat*crit[i]*sqrt(1/tau)/sqrt(n)
}
mean(upper1[1,]-lower1[1,])
mean(upper2[1,]-lower2[1,])
mean(upper3[1,]-lower3[1,])
plot(beta_hat1,type="l",ylim=c(-6,6),ylab=expression(beta[1](t)))
lines(upper1[1,],col="blue")
lines(lower1[1,],col="blue")
abline(h=0,col="red")

