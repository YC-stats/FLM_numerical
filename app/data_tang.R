rm(list=ls())
source("gcv_lamb.R")
source("loocv.R")

library(MASS)
library(ggplot2)
library(forecast)
library(abind)

time=10000
alpha=c(0.05,0.1)
kk = length(alpha)

xx1=as.matrix(read.csv("curve_demand_2015.csv", header=FALSE))
xx2=as.matrix(read.csv("curve_demand_2016.csv", header=FALSE))
xx3=as.matrix(read.csv("curve_demand_2017.csv", header=FALSE))
y1=as.matrix(read.csv("mean_demand_2015.csv", header=FALSE))
y2=as.matrix(read.csv("mean_demand_2016.csv", header=FALSE))
y3=as.matrix(read.csv("mean_demand_2017.csv", header=FALSE))
xx=cbind(xx1,xx2,xx3) #res*nn

nn=dim(xx)[2]
y_med=as.vector(rbind(y1,y2,y3))
y=log(y_med[4:nn]) #4 when using j=1:3 [4:n]
prec=as.vector(y-mean(y)) # centering

res=dim(xx)[1]

xx_1 <- xx[,-c(1,2,nn)] #X_{i1}(t) #res*(nn-3) #
xx_2 <- xx[,-c(1,nn-1,nn)]#X_{i2}(t) #
xx_3 <- xx[,-c(nn-2,nn-1,nn)]#X_{i3}(t)

temp1 <- xx_1-rowMeans(xx_1)
temp2 <- xx_2-rowMeans(xx_2)
temp3 <- xx_3-rowMeans(xx_3)
temp <- cbind(temp1,temp2,temp3) #res*(3n)
n = dim(temp1)[2]

# cov1 <- temp1%*%t(temp1)/(nn-3)
# cov2 <- temp2%*%t(temp2)/(nn-3)
# cov3 <- temp3%*%t(temp3)/(nn-3)
# 
# write.csv(cov1, "cov1.csv", row.names = FALSE)
# write.csv(cov2, "cov2.csv", row.names = FALSE)
# write.csv(cov3, "cov3.csv", row.names = FALSE)

upper=array(0,c(kk,res,3))
lower=array(0,c(kk,res,3))

eig_val_1 <- as.matrix(read.csv("eigenvalues_1.csv",header = FALSE))
eig_func_1 <- as.matrix(read.csv("eigenfunctions_1.csv", header = FALSE))
eig_val_2 <- as.matrix(read.csv("eigenvalues_2.csv",header = FALSE))
eig_func_2 <- as.matrix(read.csv("eigenfunctions_2.csv", header = FALSE))
eig_val_3 <- as.matrix(read.csv("eigenvalues_3.csv",header = FALSE))
eig_func_3 <- as.matrix(read.csv("eigenfunctions_3.csv", header = FALSE))

l <- dim(eig_func_1)[1]
count <- rep(0,2)
a_1 <- 1-1/sqrt(2)
a_2 <- 1+sqrt(2)
width = rep(0,2)
t=1:144/144
  
rho_med <- abind(eig_val_1,eig_val_2,eig_val_3)
eigen_phi_med <- array(data=c(eig_func_1,eig_func_2,eig_func_3), dim=c(dim(eig_func_1)[1],dim(eig_func_1)[2],3))

v <- 6# loocv(rho_med,eigen_phi_med,temp,prec,l,3)
rho <- c(eig_val_1[1:v,],eig_val_2[1:v,],eig_val_3[1:v,]) # 3v*1
eigen_phi <- t(rbind(eig_func_1[1:v,],eig_func_2[1:v,],eig_func_3[1:v,])) # res*3v

score <- c()  
for (i in 1:3){
  score <- cbind(score, t(temp[,((i-1)*n+1):(i*n)])%*%eigen_phi[,((i-1)*v+1):(i*v)]/res) # n*(3v)
}

D <- diag(rho) # (3v)*(3v)
  
lamb_hat <- gcv_lamb(prec,score,3*v,D)$lamb #5.07e-14

## estimation
theta_hat <- solve(t(score)%*%score/n+lamb_hat*D)%*%t(score)%*%prec/n # (3v)*1
beta_hat <- c()

for (i in 1:3){
  beta_hat <- c(beta_hat,eigen_phi[,((i-1)*v+1):(i*v)]%*%theta_hat[((i-1)*v+1):(i*v)]) #(3*res)*1
}

beta_boots <- matrix(0,3*res,time)
## bootstrap procedure
for (i in 1:time){
  boots_seq <- a_1 + (a_2 - a_1)*rbinom(n, size = 1, prob = 1/3) # n*1
  #lamb <- gcv_b(y,score,D,boots_seq)$lamb
  M <- diag(boots_seq)
  boots_theta_hat <- solve(t(score)%*%M%*%score/n+lamb_hat*D)%*%t(score)%*%M%*%prec/n
  for (k in 1:3){
    beta_boots[((k-1)*res+1):(k*res),i] <- 
      eigen_phi[,((k-1)*v+1):(k*v)]%*%boots_theta_hat[((k-1)*v+1):(k*v)]
  }
}
  
q <- sapply(1:3,function(i){
  apply(abs(sweep(beta_boots[((i-1)*res+1):(i*res),],1,beta_hat[((i-1)*res+1):(i*res)],"-")), 2, max)
  },simplify = TRUE) # time*3
ss <- apply(q,2,sort)
crit <- ss[c(floor(0.9*time), floor(0.95*time)),]

for (k in 1:3){
  for (i in 1:2){
    upper[i,,k]=beta_hat[((k-1)*res+1):(k*res)]+crit[i,k]
    lower[i,,k]=beta_hat[((k-1)*res+1):(k*res)]-crit[i,k]
  }
}

mean(2*crit[2,1])
mean(2*crit[2,2])
mean(2*crit[2,3])

data_1=data.frame(t,beta_hat=beta_hat[1:res],beta_true=rep(0,res),upper=upper[2,,1],lower=lower[2,,1])
data_2=data.frame(t,beta_hat=beta_hat[(res+1):(2*res)],beta_true=rep(0,res),upper=upper[2,,2],lower=lower[2,,2])
data_3=data.frame(t,beta_hat=beta_hat[(2*res+1):(3*res)],beta_true=rep(0,res),upper=upper[2,,3],lower=lower[2,,3])

  ggplot(data_3,aes(x=t))+
  geom_line(aes(y=beta_hat),color="black",
            linetype=1,size=1)+geom_line(aes(y=upper),
          color="blue",linetype=4,size=1)+geom_line(
          aes(y=lower),color="blue",linetype=4,size=1)+
    geom_line(aes(y=beta_true),linetype=1)+
  labs(x=expression(t),y=expression(hat(beta)[3](t)))+
  theme(legend.title=element_blank())+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15, angle=45))+
  theme(axis.title.y=element_text(vjust = 0.5,
                                  hjust = 0.5,angle =360))+
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  ggtitle("RKHS")+theme(plot.title = element_text(hjust=0.5))
#ggtitle("Mean temperature")+
#theme(plot.title = element_text(hjust = 0.5))
ggplot(data2,aes(x=t))+
  geom_line(aes(y=beta_hat2),color="black",
            linetype=1,size=1)+
  geom_line(aes(y=upper2[1,]),
            color="blue",linetype=4,size=1)+geom_line(
              aes(y=lower2[1,]),color="blue",linetype=4,size=1)+
  geom_line(aes(y=upper5[1,]),color="red",
            linetype=2,size=1)+geom_line(aes(y=lower5[1,]),
                                         color="red",linetype=2,size=1)+
  labs(x=expression(t),y=expression(hat(beta)[2](t)))+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15, angle=45))+
  geom_line(aes(y=beta_true),linetype=1)+
  theme(axis.title.y=element_text(vjust = 0.5,
                                  hjust = 0.5,angle =360))+
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  ggtitle("Legendre")+theme(plot.title = element_text(hjust=0.5))
ggplot(data3,aes(x=t))+
  geom_line(aes(y=beta_hat3),linetype=1,size=1)+
  geom_line(aes(y=upper3[1,]),
            color="blue",linetype=4,size=1)+geom_line(
              aes(y=lower3[1,]),color="blue",linetype=4,size=1)+
  geom_line(aes(y=upper6[1,]),color="red",
            linetype=2,size=1)+geom_line(aes(y=lower6[1,]),
                                         color="red",linetype=2,size=1)+
  labs(x=expression(t),y=expression(hat(beta)[3](t)))+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15, angle=45))+
  geom_line(aes(y=beta_true),linetype=1)+
  theme(axis.title.y=element_text(vjust = 0.5,
                                  hjust = 0.5,angle =360))+
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  ggtitle("Legendre")+theme(plot.title = element_text(hjust=0.5))
# ggplot(data4,aes(x=t))+
#   geom_line(aes(y=beta_hat4),linetype=1,size=1)+
#   ylim(-0.175,0.175)+geom_line(aes(y=upper4[1,]),
#   color="blue",linetype=2,size=1)+geom_line(
#   aes(y=lower4[1,]),color="blue",linetype=2,size=1)+
#   labs(x=expression(t),y=expression(beta[4](t)))
# ggplot(data5,aes(x=t))+
#   geom_line(aes(y=beta_hat5),linetype=1,size=1)+
#   ylim(-0.175,0.175)+geom_line(aes(y=upper5[1,]),
#                                color="blue",linetype=2,size=1)+geom_line(
#                                  aes(y=lower5[1,]),color="blue",linetype=2,size=1)+
#   labs(x=expression(t),y=expression(beta[5](t)))
# A=matrix(0,res,4)
# A[,1]=leg_true[1,]#cos(2*pi*t)
# A[,2]=leg_true[2,]#sin(2*pi*t)
# A[,3]=leg_true[3,]#cos(4*pi*t)
# A[,4]=leg_true[4,]#sin(4*pi*t)
# A[,5]=cos(6*pi*t)
# A[,6]=sin(6*pi*t)
# A[,7]=cos(8*pi*t)
# A[,8]=sin(8*pi*t)
# A[,9]=rep(1,res)#cos(10*pi*t)
# A[,10]=sin(10*pi*t)
# A[,11]=rep(1,res)
# result=solve(t(A)%*%A)%*%(t(A)%*%beta_hat1)
# beta_hat=A%*%result
#plot(upper1[1,],type="l",ylim=c(-0.5,0.5))
#lines(lower1[1,])
#lines(beta_hat,col="red")
# data=data.frame(t,beta_hat,upper1[1,],lower1[1,])
# ggplot(data,aes(x=t))+
#   geom_line(aes(y=beta_hat),color="black",
#             linetype=1)+ylim(-0.125,0.175)+
#   geom_line(aes(y=upper1[1,]),color="blue",linetype=2)+
#   geom_line(aes(y=lower1[1,]),color="blue",linetype=2)+
#   labs(x=expression(t),y=expression(hat(beta)[1](t)))+
#   theme(legend.title=element_blank())+
#   theme(axis.title.y=element_text(vjust = 0.5,
#                                   hjust = 0.5,angle =360))+
#   theme(axis.text.x=element_text(size=10,face="bold"))+
#   theme(axis.text.y=element_text(size=10,face="bold"))+
#   ggtitle("k=4")+theme(plot.title = element_text(hjust=0.5))




