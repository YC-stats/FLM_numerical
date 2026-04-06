rm(list=ls())
source("pbootsfpc.R")
source("mv.R")
source("gcv_lamb.R")
source("leg_constru.R")
library(MASS)
library(ggplot2)
time=10000
alpha=c(0.05,0.1)
kk=length(alpha)
xx1=as.matrix(read.csv("curve_demand_2015.csv", header=FALSE))
xx2=as.matrix(read.csv("curve_demand_2016.csv", header=FALSE))
xx3=as.matrix(read.csv("curve_demand_2017.csv", header=FALSE))
y1=as.matrix(read.csv("mean_demand_2015.csv", header=FALSE))
y2=as.matrix(read.csv("mean_demand_2016.csv", header=FALSE))
y3=as.matrix(read.csv("mean_demand_2017.csv", header=FALSE))
xx=cbind(xx1,xx2,xx3)
n=dim(xx)[2]
y_med=as.vector(rbind(y1,y2,y3))
#as.matrix(read.csv("price.csv", header=FALSE))[730:(730+n),]
xx=xx[,-n]
y=log(y_med[4:n])
#log(as.numeric(y_med[4:n,7]))
nn=dim(xx)[2]
res=dim(xx)[1]
xx_1=xx[,-c(1,2)] #X_{i1}(t) #res*n1
xx_2=xx[,-c(1,nn)]#X_{i2}(t)
xx_3=xx[,-c((nn-1),nn)]#X_{i3}(t)
x_med1=prcomp(t(xx_1),rank=10,center=TRUE)
x_med2=prcomp(t(xx_2),rank=10,center=TRUE)
x_med3=prcomp(t(xx_3),rank=10,center=TRUE)
cvp=0
mm=0
while (cvp <= 0.95){
  mm=mm+1
  cvp=summary(x_med1)$importance[3*mm]
}
c1=2*mm#floor(sqrt(log(dim(xx_1)[2]))*mm)
cvp=0
mm=0
while (cvp <= 0.95){
  mm=mm+1
  cvp=summary(x_med2)$importance[3*mm]
}
c2=2*mm#floor(sqrt(log(dim(xx_2)[2]))*mm)
cvp=0
mm=0
while (cvp <= 0.95){
  mm=mm+1
  cvp=summary(x_med3)$importance[3*mm]
}
c3=2*mm#floor(sqrt(log(dim(xx_3)[2]))*mm)
c=c1+c2+c3
n=dim(xx_1)[2]
t=(1:res)/res
x1=matrix(0,n,c1)
x2=matrix(0,n,c2)
x3=matrix(0,n,c3)
temp1=matrix(0,res,n)
temp2=matrix(0,res,n)
temp3=matrix(0,res,n)
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
mu=mean(y)
mu_1=apply(xx_1,1,mean) #res*1 vector
mu_2=apply(xx_2,1,mean)
mu_3=apply(xx_3,1,mean)
for (i in 1:res){#centering
  temp1[i,]=xx_1[i,]-mu_1[i]#res*n matrix
  temp2[i,]=xx_2[i,]-mu_2[i]
  temp3[i,]=xx_3[i,]-mu_3[i]
}
c_max=max(c1,c2,c3)
basis=matrix(1,c_max,res)
basis[2:c_max,]=sqrt(2)*cos(pi*(1:(c_max-1))%*%t(t))
coef_1=t(ginv(t(basis))%*%temp1) #n*c1
coef_2=t(ginv(t(basis))%*%temp2)
coef_3=t(ginv(t(basis))%*%temp3)
deri_2=c(0,-sqrt(2)*(pi*(1:(c_max-1)))^2)
deri_x1=coef_1%*%sweep(basis,1,deri_2,"*") #n*res
deri_x2=coef_2%*%sweep(basis,1,deri_2,"*")
deri_x3=coef_3%*%sweep(basis,1,deri_2,"*")
prec=as.vector(y-mu)#centering
aa1=prcomp(t(temp1),rank=c1)
aa2=prcomp(t(temp2),rank=c2)
aa3=prcomp(t(temp3),rank=c3)
pcs1=aa1$x/sqrt(res)  #n*c1 matrix 
pcs2=aa2$x/sqrt(res) 
pcs3=aa3$x/sqrt(res) 
eig1=as.matrix(aa1$rotation*sqrt(res)) 
eig2=as.matrix(aa2$rotation*sqrt(res)) 
eig3=as.matrix(aa3$rotation*sqrt(res))
x_deri_fpc1=t(ginv(pcs1)%*%deri_x1)# res*c_max
x_deri_fpc2=t(ginv(pcs2)%*%deri_x2)
x_deri_fpc3=t(ginv(pcs3)%*%deri_x3)
f1=apply(pcs1,2,sd) #a dimension
f2=apply(pcs2,2,sd)
f3=apply(pcs3,2,sd)
x_scale1=sweep(pcs1,2,f1,"/")#n*c1 matrix
x_scale2=sweep(pcs2,2,f2,"/")
x_scale3=sweep(pcs3,2,f3,"/")
x=cbind(x_scale1,x_scale2,x_scale3)#n*c matrix
P_med1=t(x_deri_fpc1[,1:c1])%*%x_deri_fpc1[,1:c1]/res
P_med2=t(x_deri_fpc2[,1:c2])%*%x_deri_fpc2[,1:c2]/res
P_med3=t(x_deri_fpc3[,1:c3])%*%x_deri_fpc3[,1:c3]/res
P1=sweep(sweep(P_med1,1,f1[1:c1],"/"),2,f1[1:c1],"/")
P2=sweep(sweep(P_med2,1,f2[1:c2],"/"),2,f2[1:c2],"/")
P3=sweep(sweep(P_med3,1,f3[1:c3],"/"),2,f3[1:c3],"/")
R_med=matrix(0,c,c)
R_med[1:c1,1:c1]=P1
R_med[(c1+1):(c1+c2),(c1+1):(c1+c2)]=P2
R_med[(c1+c2+1):c,(c1+c2+1):c]=P3
lamb=gcv_lamb(prec,x,c,R_med)$lamb
R=lamb*R_med
theta_hat=as.vector(solve(t(x)%*%x/n+R)%*%t(x)%*%prec/n)
f=c(f1,f2,f3)#(5a)*1 vector
b_hat=theta_hat/f#(5a)*1 vector
beta_hat1=as.vector(eig1%*%b_hat[1:c1])#res*1 vector
beta_hat2=as.vector(eig2%*%b_hat[(c1+1):(c1+c2)])
beta_hat3=as.vector(eig3%*%b_hat[(c1+c2+1):(c1+c2+c3)])
eps=prec-x%*%theta_hat
for (i in 1:n)
  z[,i]=t(x[i,]*eps[i])
S=t(apply(z,1,cumsum))##(5a)*n matrix
beta_true=rep(0,res)
####choose window size
m=mv(S,range_m=(floor(dim(S)[2]^(1/3)*0.25)):
       (floor(2*dim(S)[2]^(1/3))))$m
####construct SCB#####
#(S,x,list(f1,f2,f3),t,m,c(c1,c2,c3),time,R,leg_true)
med=pbootsfpc(S,x,list(f1,f2,f3),t,m,c(c1,c2,c3),time,
              R,list(eig1,eig2,eig3))
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
#data4=data.frame(t,beta_hat4,upper4[1,],lower4[1,])
#data5=data.frame(t,beta_hat5,upper5[1,],lower5[1,])
#g1 ylim=c(-0.1,0.1)
ggplot(data1,aes(x=t))+
  geom_line(aes(y=beta_hat1),color="black",linetype=1,size=1)+
  geom_line(aes(y=upper1[1,]),
  color="blue",linetype=4,size=1)+geom_line(
  aes(y=lower1[1,]),color="blue",linetype=4,size=1)+
  geom_line(aes(y=upper4[1,]),color="red",
  linetype=2,size=1)+geom_line(aes(y=lower4[1,]),
  color="red",linetype=2,size=1)+
  labs(x=expression(t),y=expression(hat(beta)[1](t)))+
  theme(legend.title=element_blank())+
  geom_line(aes(y=beta_true),linetype=1)+
  theme(axis.title.y=element_text(vjust = 0.5,
      hjust = 0.5,angle =360,size=15),
      axis.title.x=element_text(size=15))+
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  ggtitle("FPC")+theme(plot.title = element_text(hjust=0.5))
#ggtitle("Mean temperature")+
#theme(plot.title = element_text(hjust = 0.5))
ggplot(data2,aes(x=t))+
  geom_line(aes(y=beta_hat2),linetype=1,size=1)+
  geom_line(aes(y=upper2[1,]),
  color="blue",linetype=4,size=1)+geom_line(
  aes(y=lower2[1,]),color="blue",linetype=4,size=1)+
  geom_line(aes(y=upper5[1,]),color="red",
  linetype=2,size=1)+geom_line(aes(y=lower5[1,]),
  color="red",linetype=2,size=1)+
  labs(x=expression(t),y=expression(hat(beta)[2](t)))+
  theme(legend.title=element_blank())+
  geom_line(aes(y=beta_true),linetype=1)+
  theme(axis.title.y=element_text(size=15,vjust = 0.5,
  hjust = 0.5,angle =360),
  axis.title.x=element_text(size=15))+
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  ggtitle("FPC")+theme(plot.title = element_text(hjust=0.5))
ggplot(data3,aes(x=t))+ylim(-8,10)+
  geom_line(aes(y=beta_hat3),linetype=1,size=1)+
  geom_line(aes(y=upper3[1,]),
  color="blue",linetype=4,size=1)+geom_line(
  aes(y=lower3[1,]),color="blue",linetype=4,size=1)+
  geom_line(aes(y=upper6[1,]),color="red",
  linetype=2,size=1)+geom_line(aes(y=lower6[1,]),
  color="red",linetype=2,size=1)+
  labs(x=expression(t),y=expression(hat(beta)[3](t)))+
  theme(legend.title=element_blank())+
  geom_line(aes(y=beta_true),linetype=1)+
  theme(axis.title.y=element_text(vjust = 0.5,
    hjust = 0.5,angle =360))+
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  ggtitle("FPC")+theme(plot.title = element_text(hjust=0.5))




