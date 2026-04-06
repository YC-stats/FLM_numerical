rm(list=ls())
source("pbootslag.R")
source("mv.R")
source("gcv_lamb.R")
source("leg_constru.R")
library(MASS)
library(ggplot2)
library(forecast)

time=10000
set.seed(123)
alpha=c(0.05,0.1)
kk=length(alpha)

## Upload datasets ###
xx1=as.matrix(read.csv("curve_demand_2015.csv", header=FALSE))
xx2=as.matrix(read.csv("curve_demand_2016.csv", header=FALSE))
xx3=as.matrix(read.csv("curve_demand_2017.csv", header=FALSE))
y1=as.matrix(read.csv("mean_demand_2015.csv", header=FALSE))
y2=as.matrix(read.csv("mean_demand_2016.csv", header=FALSE))
y3=as.matrix(read.csv("mean_demand_2017.csv", header=FALSE))
temp_z= read.csv("daily_temp.csv") # NA

temp_z <- temp_z$national_tmax

## select the cutoff c by searching over possible thresholds and choosing the one that 
## maximizes the absolute correlation between the response and indicator predictor

basis_name <- "Legen"

xx=cbind(xx1,xx2,xx3) # res*nn
y_med=as.vector(rbind(y1,y2,y3)) # nn*1
nn=dim(xx)[2] # column/sample size
res=dim(xx)[1] # row 

## centering
y=log(y_med[3:nn]) #4 when using j=1:3 [4:n]
prec=as.vector(y-mean(y)) 

## centered predictors
l <- 2
#pred_curve <- array(data = c(xx[,-c(1,2,nn)], xx[,-c(1,nn-1,nn)], xx[,-c(nn-2,nn-1,nn)]), dim = c(res,nn-l,l))
pred_curve <- array(data = c(xx[,-c(1,nn)], xx[,-c(nn-1,nn)]), dim = c(res,nn-l,l))
row_means <- apply(pred_curve, c(1,3), mean)
x <- sweep(pred_curve, c(1,3), row_means, "-")

if (!any(is.na(temp_z))){
  zz <- temp_z[-c(1,2)]
  z_temp_1 <- as.numeric(zz >= 31.7)
  z_temp_2 <- as.numeric(zz <= 14.6) # c()
  z_temp <- cbind(z_temp_1,z_temp_2)
}

temp_dim <- dim(z_temp)[2]

## candidate cutoffs
# cand <- sort(unique(zz))
# 
# cors <- sapply(cand, function(c){
#   z <- as.numeric(zz >= c)
#   if (sd(z) == 0) return(0)
#   cor(y, z)
# })
# 
# c_star <- cand[which.max(cors)] # or which.min
## index for heating: 148 from 780 series; index for cooling: 298

## truncation number
c_n <- rep(0,l)
for (i in 1:l){
  x_med=prcomp(t(x[,,i]),rank=10,center=TRUE)
  cvp=0
  mm=0
  while (cvp <= 0.95){
    mm=mm+1
    cvp=summary(x_med)$importance[3*mm]
  }
  c_n[i]=2*mm
}

c=sum(c_n)
c_max <- max(c_n)
n=dim(x[,,1])[2]
t=(1:res)/res

eps=rep(0,n)

beta_hat <- matrix(0,res,l)

#######fit model#######
score <- array(0,dim=c(n,c_n[1],l))
x_scale <- array(0,dim=c(n,c_n[1],l))
f <- matrix(0,c_n[1],l)
R_med=matrix(0,c,c)

## Legendre basis ##
if (basis_name == "Legen"){
  leg_true=matrix(0,nrow=c_max,ncol=res)
  tilde_t=2*t-1
  for (i in 1:c_max)
    leg_true[i,]=sqrt(2*i-1)*leg_constru(tilde_t,i-1)
  
  for (i in 1:l){
    score[,,i] <- t(ginv(t(leg_true[1:c_n[i],]))%*%x[,,i]) #n*c
    f[,i] <- apply(score[,,i],2,sd) #c dimension
    x_scale[,,i] <- sweep(score[,,i],2,f[,i],"/")
  }
  
  P=diag(c(0,0,720,8400,49680,203280)/(f[,1]^2),c_n[1],c_n[1])
  P[3,5]=1440*sqrt(5)/f[3,1]/f[5,1]
  P[5,3]=P[3,5]
  P[4,6]=3360*sqrt(77)/f[4,1]/f[6,1]
  P[6,4]=P[4,6]
  
  for (k in 1:l){
    R_med[((k-1)*c_n[k]+1):(k*c_n[k]),((k-1)*c_n[k]+1):(k*c_n[k])]=P
  }
}else{
  ## FPC basis ##
  basis_spline <- fda::create.bspline.basis(rangeval=c(0,1), nbasis=25, norder=4)
  fpc_basis <- array(0,c(res,c_n[1],l))
  harmonics_2nd <- array(0,c(res,c_n[1],l))
  for (i in 1:l){
    fd.data <- fda::smooth.basis(1:res/res, x[,,i], basis_spline)$fd
    fd_pca <- fda::pca.fd(fd.data, nharm=c_n[i], fda::fdPar(fd.data), centerfns=TRUE)
    score[,,i] <- fd_pca$scores # n*c
    f[,i] <- apply(score[,,i],2,sd) 
    x_scale[,,i] <- sweep(score[,,i],2,f[,i],"/")
    fpc_basis[,,i] <- fda::eval.fd(1:res/res, fd_pca$harmonics)
    harmonics_2nd[,,i] <- fda::eval.fd(1:res/res,fda::deriv.fd(fd_pca$harmonics, deriv = 2))
  }
  for (k in 1:l){
    R_med[((k-1)*c_n[k]+1):(k*c_n[k]),((k-1)*c_n[k]+1):(k*c_n[k])] <- 
      sweep(sweep(t(harmonics_2nd[,,k])%*%harmonics_2nd[,,k]/res,1,f[,k],"/"),2,f[,k],"/")
  }
}

x_score <- matrix(x_scale, nrow = n)

if (any(is.na(temp_z))){
  lamb=gcv_lamb(prec,x_score,c,R_med)$lamb
  R=lamb*R_med
  theta_hat=as.vector(solve(t(x_score)%*%x_score/n+R)%*%t(x_score)%*%prec/n)
  b_hat=theta_hat/as.vector(f) 
  eps=as.vector(prec-x_score%*%theta_hat)
  z=matrix(0,c,n)
  crit <- matrix(0,kk,l)
  upper <- array(0, c(kk,res,l))
  lower <- array(0, c(kk,res,l))
  upper_p <- array(0, c(kk,res,l))
  lower_p <- array(0, c(kk,res,l))
}else{
  crit <- matrix(0,kk,l+temp_dim)
  upper <- array(0, c(kk,res,l+temp_dim))
  lower <- array(0, c(kk,res,l+temp_dim))
  upper_p <- array(0, c(kk,res,l+temp_dim))
  lower_p <- array(0, c(kk,res,l+temp_dim))
  z=matrix(0,c+temp_dim,n)
  x_score <- cbind(x_score,z_temp)
  R_med <- cbind(rbind(R_med,matrix(0,temp_dim,c)),matrix(0,c+temp_dim,temp_dim))
  lamb=gcv_lamb(prec,x_score,c+temp_dim,R_med)$lamb # needs to check
  R=lamb*R_med
  theta_hat <- as.vector(solve(t(x_score)%*%x_score/n+R)%*%t(x_score)%*%prec/n)
  b_hat <- theta_hat[-c((c+1):(c+temp_dim))]/as.vector(f) 
  b_cov <- theta_hat[-c(1:c)]
  eps=as.vector(prec-x_score%*%theta_hat)
}

if (basis_name == "Legen"){
  for (i in 1:l)
    beta_hat[,i] <- as.vector(t(leg_true[1:c_n[i],])%*%b_hat[((i-1)*c_n[i]+1):(i*c_n[i])])
  leg =  rep(list(leg_true), l)
}else{
  for (i in 1:l)
    beta_hat[,i] <- as.vector(fpc_basis[,,i]%*%b_hat[((i-1)*c_n[i]+1):(i*c_n[i])])
  leg = lapply(1:l, function(i) t(fpc_basis[,,i]))
}


for (i in 1:n)
  z[,i]=x_score[i,]*eps[i]

S=t(apply(z,1,cumsum)) 

### choose window size ###
m=mv(S,range_m=(floor(dim(S)[2]^(1/3)*0.25)):
       (floor(2*dim(S)[2]^(1/3))))$m

### construct SCB ###
med=pbootslag(S,x_score,f,t,m,a=c_n,
              time,R,leg,z=z_temp) #  z=z_temp

r <- apply(med$distr,2,sort)

for (i in 1:dim(r)[2]){
  crit[,i] <- r[floor(time*(c(1,1)-alpha)),i]
}

for (i in 1:l){
  upper[,,i]=sweep(crit[,i]%*%t(med$g$g_hat[,i])/sqrt(n),2,beta_hat[,i],"+")
  lower[,,i]=sweep(-crit[,i]%*%t(med$g$g_hat[,i])/sqrt(n),2,beta_hat[,i],"+")
  upper_p[,,i]=sweep(c(1.96, 1.65)%*%t(med$g_or[,i])/sqrt(n),2,beta_hat[,i],"+")
  lower_p[,,i]=sweep(-c(1.96, 1.65)%*%t(med$g_or[,i])/sqrt(n),2,beta_hat[,i],"+")
}

upper_cov=sweep(sweep(crit[,(l+1):(l+temp_dim),drop=FALSE],2,med$g$g_theta/sqrt(n),"*"), 2, b_cov, "+")
lower_cov=sweep(sweep(-crit[,(l+1):(l+temp_dim),drop=FALSE],2,med$g$g_theta/sqrt(n),"*"), 2, b_cov, "+")
#upper_p[,,(l+1):(l+temp_dim)]=c(1.96, 1.65)*med$g$g_theta/sqrt(n)+b_cov
#lower_p[,,(l+1):(l+temp_dim)]=-c(1.96, 1.65)*med$g$g_theta/sqrt(n)+b_cov

data1=data.frame(t,beta_hat=beta_hat[,1],beta_true=rep(0,res),
                 upper=upper[1,,1],lower=lower[1,,1],
                 upper_p=upper_p[1,,1],lower_p=lower_p[1,,1])
data2=data.frame(t,beta_hat=beta_hat[,2],beta_true=rep(0,res),
                 upper=upper[1,,2],lower=lower[1,,2],
                 upper_p=upper_p[1,,2],lower_p=lower_p[1,,2])
data3=data.frame(t,beta_hat=beta_hat[,3],beta_true=rep(0,res),
                 upper=upper[1,,3],lower=lower[1,,3],
                 upper_p=upper_p[1,,3],lower_p=lower_p[1,,3])

ggplot(data2,aes(x=t))+
  geom_line(aes(y=beta_hat),color="black",
            linetype=1,size=1)+geom_line(aes(y=upper),
            color="blue",linetype=4,size=1)+geom_line(
          aes(y=lower),color="blue",linetype=4,size=1)+
  geom_line(aes(y=upper_p),color="red",
          linetype=2,size=1)+geom_line(aes(y=lower_p),
          color="red",linetype=2,size=1)+
  labs(x=expression(t),y=expression(hat(beta)[1](t)))+
  theme(legend.title=element_blank())+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15, angle=45))+
  geom_line(aes(y=beta_true),linetype=1)+
  theme(axis.title.y=element_text(vjust = 0.5,
                                  hjust = 0.5,angle =360))+
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  ggtitle("Legendre")+theme(plot.title = element_text(hjust=0.5))
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
