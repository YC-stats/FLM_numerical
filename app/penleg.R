rm(list=ls())
source("pbootslag.R")
source("mv.R")
source("gcv_lamb.R")
source("leg_constru.R")
library(MASS)
library(ggplot2)
library(forecast)
time=10000
alpha=c(0.05,0.1)
kk=length(alpha)
xx1=as.matrix(read.csv("curve_demand_2015.csv", header=FALSE))
xx2=as.matrix(read.csv("curve_demand_2016.csv", header=FALSE))
xx3=as.matrix(read.csv("curve_demand_2017.csv", header=FALSE))
y1=as.matrix(read.csv("mean_demand_2015.csv", header=FALSE))
y2=as.matrix(read.csv("mean_demand_2016.csv", header=FALSE))
y3=as.matrix(read.csv("mean_demand_2017.csv", header=FALSE))
xx=cbind(xx1,xx2,xx3) #res*nn
nn=dim(xx)[2]
y_med=as.vector(rbind(y1,y2,y3))
  #as.matrix(read.csv("price.csv", header=FALSE))[730:(730+n),]
#xx=xx[,-n] #remove the last item
y=log(y_med[4:nn]) #4 when using j=1:3 [4:n]
prec=as.vector(y-mean(y)) # centering

### arima model prediction ###
# mse <- rep(0,7)
 # for (i in (n_total-7):(n_total-1)){
 #   model <- auto.arima(prec[1:i], max.p=8, ic="aic", allowmean = FALSE) #prec[1:i]
#   mse[i-n_total+8] <- (as.numeric(forecast(model, h = 1)$mean) - prec[i+1])^2
# }
# mean(mse) #0.01373204

res=dim(xx)[1]

xx_1 <- xx[,-c(1,2,nn)] #X_{i1}(t) #res*n1 #
xx_2 <- xx[,-c(1,nn-1,nn)]#X_{i2}(t) #
xx_3=xx[,-c(nn-2,nn-1,nn)]#X_{i3}(t)

  
x_med1=prcomp(t(xx_1),rank=10,center=TRUE)
x_med2=prcomp(t(xx_2),rank=10,center=TRUE)
x_med3=prcomp(t(xx_3),rank=10,center=TRUE)
cvp=0
mm=0
while (cvp <= 0.95){
  mm=mm+1
  cvp=summary(x_med1)$importance[3*mm]
}
c1=2*mm #floor(sqrt(log(dim(xx_1)[2]))*mm)
cvp=0
mm=0
while (cvp <= 0.95){
  mm=mm+1
  cvp=summary(x_med2)$importance[3*mm]
  }
c2=2*mm
cvp=0
mm=0
while (cvp <= 0.95){
  mm=mm+1
  cvp=summary(x_med3)$importance[3*mm]
  }
c3=2*mm
c=c1+c2+c3
n=dim(xx_1)[2]
t=(1:res)/res
x1=matrix(0,n,c1)
x2=matrix(0,n,c2)
x3=matrix(0,n,c3)

eps=rep(0,n)
z=matrix(0,c,n)
upper1=matrix(0,kk,res)
lower1=matrix(0,kk,res)
upper2=matrix(0,kk,res)
lower2=matrix(0,kk,res)
upper3=matrix(0,kk,res)
lower3=matrix(0,kk,res)
upper4=matrix(0,kk,res)
lower4=matrix(0,kk,res)
upper5=matrix(0,kk,res)
lower5=matrix(0,kk,res)
upper6=matrix(0,kk,res)
lower6=matrix(0,kk,res)

#######fit model#######
temp1 <- xx_1-rowMeans(xx_1)
temp2 <- xx_2-rowMeans(xx_2)
temp3 <- xx_3-rowMeans(xx_3)

c_max=max(c1,c2,c3)
leg_true=matrix(0,nrow=c_max,ncol=res)
tilde_t=2*t-1
for (i in 1:c_max)
  leg_true[i,]=sqrt(2*i-1)*leg_constru(tilde_t,i-1)

x1=ginv(t(leg_true[1:c1,]))%*%temp1 #c*n
x2=ginv(t(leg_true[1:c2,]))%*%temp2
x3=ginv(t(leg_true[1:c3,]))%*%temp3

f1=apply(x1,1,sd) #c dimension
f2=apply(x2,1,sd)
f3=apply(x3,1,sd)
x_scale1=sweep(t(x1),2,f1,"/") #n*c matrix
x_scale2=sweep(t(x2),2,f2,"/")
x_scale3=sweep(t(x3),2,f3,"/")
x=cbind(x_scale1,x_scale2,x_scale3) #n*c matrix
P1=diag(c(0,0,720,8400,49680,203280)/(f1^2),c1,c1)
P1[3,5]=1440*sqrt(5)/f1[3]/f1[5]
P1[5,3]=1440*sqrt(5)/f1[3]/f1[5]
P1[4,6]=3360*sqrt(77)/f1[4]/f1[6]
P1[6,4]=3360*sqrt(77)/f1[4]/f1[6]
P2=diag(c(0,0,720,8400,49680,203280)/(f2^2),c2,c2)
P2[3,5]=1440*sqrt(5)/f2[3]/f2[5]
P2[5,3]=1440*sqrt(5)/f2[3]/f2[5]
P2[4,6]=3360*sqrt(77)/f2[4]/f2[6]
P2[6,4]=3360*sqrt(77)/f2[4]/f2[6]
P3=diag(c(0,0,720,8400,49680,203280)/(f3^2),c3,c3)
P3[3,5]=1440*sqrt(5)/f3[3]/f3[5]
P3[5,3]=1440*sqrt(5)/f3[3]/f3[5]
P3[4,6]=3360*sqrt(77)/f3[4]/f3[6]
P3[6,4]=3360*sqrt(77)/f3[4]/f3[6]

R_med=matrix(0,c,c)
R_med[1:c1,1:c1]=P1
R_med[(c1+1):(c1+c2),(c1+1):(c1+c2)]=P2
R_med[(c1+c2+1):c,(c1+c2+1):c]=P3

lamb=gcv_lamb(prec,x,c,R_med)$lamb
R=lamb*R_med
theta_hat=as.vector(solve(t(x)%*%x/n+R)%*%t(x)%*%prec/n)
f=c(f1,f2,f3) #(3c)*1 vector
b_hat=theta_hat/f #(3c)*1 vector
beta_hat1=as.vector(t(leg_true[1:c1,])%*%b_hat[1:c1])#res*1 vector
beta_hat2=as.vector(t(leg_true[1:c2,])%*%b_hat[(c1+1):(c1+c2)])
beta_hat3=as.vector(t(leg_true[1:c3,])%*%b_hat[(c1+c2+1):c])

eps=as.vector(prec-x%*%theta_hat)
for (i in 1:n)
  z[,i]=x[i,]*eps[i]
S=t(apply(z,1,cumsum)) ##(5a)*n matrix
beta_true=rep(0,res)
####choose window size
m=mv(S,range_m=(floor(dim(S)[2]^(1/3)*0.25)):
       (floor(2*dim(S)[2]^(1/3))))$m
####construct SCB#####
med=pbootslag(S,x,list(f1,f2,f3),t,m,a=c(c1,c2,c3),
              time,R,leg_true)
r=sort(med$distr)
g1=med$g1;g2=med$g2;g3=med$g3
g1_or=med$g1_or;g2_or=med$g2_or;g3_or=med$g3_or
crit=r[floor(time*(c(1,1)-alpha))]
for (i in 1:kk){
  upper1[i,]=beta_hat1+crit[i]*g1/sqrt(n)#crit[i]
  lower1[i,]=beta_hat1-crit[i]*g1/sqrt(n)
  upper2[i,]=beta_hat2+crit[i]*g2/sqrt(n)
  lower2[i,]=beta_hat2-crit[i]*g2/sqrt(n)
  upper3[i,]=beta_hat3+crit[i]*g3/sqrt(n)
  lower3[i,]=beta_hat3-crit[i]*g3/sqrt(n)
  upper4[i,]=beta_hat1+1.96*g1_or/sqrt(n)
  lower4[i,]=beta_hat1-1.96*g1_or/sqrt(n)
  upper5[i,]=beta_hat2+1.96*g2_or/sqrt(n)
  lower5[i,]=beta_hat2-1.96*g2_or/sqrt(n)
  upper6[i,]=beta_hat3+1.96*g3_or/sqrt(n)
  lower6[i,]=beta_hat3-1.96*g3_or/sqrt(n)
}
data1=data.frame(t,beta_hat1,upper1[1,],lower1[1,],
                 beta_true,upper4[1,],lower4[1,])#95%
data2=data.frame(t,beta_hat2,upper2[1,],lower2[1,],
                 beta_true,upper5[1,],lower5[1,])
data3=data.frame(t,beta_hat3,upper3[1,],lower3[1,],
                 beta_true,upper6[1,],lower6[1,])
ggplot(data1,aes(x=t))+
  geom_line(aes(y=beta_hat1),color="black",
      linetype=1,size=1)+geom_line(aes(y=upper1[1,]),
      color="blue",linetype=4,size=1)+geom_line(
  aes(y=lower1[1,]),color="blue",linetype=4,size=1)+
  geom_line(aes(y=upper4[1,]),color="red",
  linetype=2,size=1)+geom_line(aes(y=lower4[1,]),
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




