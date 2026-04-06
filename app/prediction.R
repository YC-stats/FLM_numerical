source("leg_constru.R")
source("gcv_lamb.R")
library(MASS)

xx1=as.matrix(read.csv("curve_demand_2015.csv", header=FALSE))
xx2=as.matrix(read.csv("curve_demand_2016.csv", header=FALSE))
xx3=as.matrix(read.csv("curve_demand_2017.csv", header=FALSE))
y1=as.matrix(read.csv("mean_demand_2015.csv", header=FALSE))
y2=as.matrix(read.csv("mean_demand_2016.csv", header=FALSE))
y3=as.matrix(read.csv("mean_demand_2017.csv", header=FALSE))
temp_z=log(read.csv("daily_temp.csv")$national_tmed)
xx=cbind(xx1,xx2,xx3)
n_total=dim(xx)[2]
y_med=as.vector(rbind(y1,y2,y3))
y=log(y_med) #3:n
mse_func <- rep(0,15) # changed to 15 for prediction interval
res=dim(xx)[1]
pred_interval <- c()
response_true <- c()
basis_name <- "Legen"

# if (!any(is.na(temp_z))){
#   zz <- temp_z$national_tmax[-c(1,nn)] # tmin or tmax
#   z_temp <- zz - mean(zz)
# }

for (k in (n_total-14):n_total){
  z_cut <- temp_z[1:k]
  x_cut <- xx[,1:k]
  nn=dim(x_cut)[2]
  xx_1=x_cut[,-c(1,nn)]
  ##Y_3= beta_1x_2+beta_2x_1
  ##Y_4=beta_1x_3+beta_2x_2
  #Y_n+1 = beta_1x_n+beta_2x_n-1
  xx_2=x_cut[,-c(nn-1,nn)]
  prec <- y[1:k]
  prec_cut <- prec[-c(1,2)]-mean(prec[-c(1,2)]) 
  
  x_med1=prcomp(t(xx_1),rank=10,center=TRUE)
  x_med2=prcomp(t(xx_2),rank=10,center=TRUE)
  cvp=0
  mm=0
  while (cvp <= 0.95){
    mm=mm+1
    cvp=summary(x_med1)$importance[3*mm]
  }
  c1=2*mm 
  cvp=0
  mm=0
  while (cvp <= 0.95){
    mm=mm+1
    cvp=summary(x_med2)$importance[3*mm]
  }
  c2=2*mm
  c_n = c(c1,c2)
  n=dim(xx_1)[2]
  t=(1:res)/res
  x1=matrix(0,n,c1)
  x2=matrix(0,n,c2)
  temp1=matrix(0,res,n)
  temp2=matrix(0,res,n)
  mu_1=apply(xx_1,1,mean) #res*1 vector
  mu_2=apply(xx_2,1,mean)
  for (i in 1:res){#centering
    temp1[i,]=xx_1[i,]-mu_1[i]#res*n matrix
    temp2[i,]=xx_2[i,]-mu_2[i]
  }
  x <- array(data = c(temp1,temp2),dim=c(res,n,2))
  z_temp <- z_cut[-c(1,nn)] - mean(z_cut[-c(1,nn)])
  
  c = c1+c2
  c_max=max(c1,c2) #,c3)
  if (basis_name == "Legen"){
    leg_true=matrix(0,nrow=c_max,ncol=res)
    tilde_t=2*t-1
    for (i in 1:c_max)
      leg_true[i,]=sqrt(2*i-1)*leg_constru(tilde_t,i-1)
    for (i in 1:n){
      x1[i,]=ginv(t(leg_true[1:c1,]))%*%temp1[,i]#n*a matrix
      x2[i,]=ginv(t(leg_true[1:c2,]))%*%temp2[,i]
    }
    
    f1=apply(x1,2,sd)#a dimension
    f2=apply(x2,2,sd)
    x_scale1=sweep(x1,2,f1,"/")#n*c matrix
    x_scale2=sweep(x2,2,f2,"/")
    x=cbind(x_scale1,x_scale2)
    
    P1=diag(c(0,0,720,8400,49680,203280)/(f1^2),c1,c1)#
    P1[3,5]=1440*sqrt(5)/f1[3]/f1[5]
    P1[5,3]=1440*sqrt(5)/f1[3]/f1[5]
    P1[4,6]=3360*sqrt(77)/f1[4]/f1[6]
    P1[6,4]=3360*sqrt(77)/f1[4]/f1[6]
    P2=diag(c(0,0,720,8400,49680,203280)/(f2^2),c2,c2)#
    P2[3,5]=1440*sqrt(5)/f2[3]/f2[5]
    P2[5,3]=1440*sqrt(5)/f2[3]/f2[5]
    P2[4,6]=3360*sqrt(77)/f2[4]/f2[6]
    P2[6,4]=3360*sqrt(77)/f2[4]/f2[6]
    R_med=matrix(0,c,c)
    R_med[1:c1,1:c1]=P1
    R_med[(c1+1):(c1+c2),(c1+1):(c1+c2)]=P2
  }else{
    R_med=matrix(0,c,c)
    score <- array(0,c(n,c_max,2))
    f <- matrix(0,c_max,2)
    x_scale <- array(0,c(n,c_max,2))
    basis_spline <- fda::create.bspline.basis(rangeval=c(0,1), nbasis=25, norder=4)
    fpc_basis <- array(0,c(res,c_max,2))
    harmonics_2nd <- array(0,c(res,c_max,2))
    for (i in 1:2){
      fd.data <- fda::smooth.basis(1:res/res, x[,,i], basis_spline)$fd
      fd_pca <- fda::pca.fd(fd.data, nharm=c_n[i], fda::fdPar(fd.data), centerfns=TRUE)
      score[,,i] <- fd_pca$scores # n*c
      f[,i] <- apply(score[,,i],2,sd) 
      x_scale[,,i] <- sweep(score[,,i],2,f[,i],"/")
      fpc_basis[,,i] <- fda::eval.fd(1:res/res, fd_pca$harmonics)
      harmonics_2nd[,,i] <- fda::eval.fd(1:res/res,fda::deriv.fd(fd_pca$harmonics, deriv = 2))
    }
    for (k in 1:2){
      R_med[((k-1)*c_n[k]+1):(k*c_n[k]),((k-1)*c_n[k]+1):(k*c_n[k])] <- 
        sweep(sweep(t(harmonics_2nd[,,k])%*%harmonics_2nd[,,k]/res,1,f[,k],"/"),2,f[,k],"/")
    }
    x <- matrix(x_scale, nrow = n)
  }
  
  lamb=gcv_lamb(prec_cut,x,c,R_med)$lamb
  R=lamb*R_med
  theta_hat=as.vector(solve(t(x)%*%x/n+R)%*%t(x)%*%prec_cut/n)
    # med_term <- solve(t(x)%*%x/n+R)
    # med_num <- solve(t(z_temp)%*%z_temp)
    # beta_cov <- as.numeric(med_num%*%(t(prec_cut)%*%z_temp-t(z_temp)%*%x%*%med_term%*%t(x)%*%prec_cut/n))/
    #  (1-as.numeric(med_num%*%t(z_temp)%*%x%*%med_term%*%t(x)%*%z_temp/n))
    # theta_hat <- as.vector(med_term%*%(t(x)%*%prec_cut/n-beta_cov*t(x)%*%z_temp/n))
  residual <- as.vector(prec_cut - x%*%theta_hat)
  res_quan <- sort(residual)[c(floor(0.005*n),floor(0.995*n))]
  pred_interval <- rbind(pred_interval,as.numeric(x[n,,drop=FALSE]%*%theta_hat)+c(res_quan[1],res_quan[2]))
  mse_func[k-n_total+15] <- (x[n,,drop=FALSE]%*%theta_hat-prec_cut[n])^2
  response_true <- c(response_true, prec_cut[length(prec_cut)])
}

mean(mse_func) # 0.00845577 decrease a little when adding temp 0.0082
library(ggplot2)
data <- data.frame(Day=c(1:15),response_true=response_true, upper=pred_interval[,1],
                   lower=pred_interval[,2])
ggplot(data = data, aes(x = Day)) + 
  geom_line(aes(y=response_true),color="black",linetype=1,size=1)+
  geom_line(aes(y=upper),color="blue",linetype=2,size=1)+
  geom_line(aes(y=lower),color="blue",linetype=2,size=1)+
  labs(x="Day",y="log(Demand)")+ylim(-0.35,0.25) + 
  theme(legend.title=element_blank())+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15, angle=45))+
  theme(axis.title.y=element_text(vjust = 0.8,
                                         hjust = 0.5,angle =90))+
         theme(axis.text.x=element_text(size=15,face="bold"))+
         theme(axis.text.y=element_text(size=15,face="bold"))+
         #ggtitle("99% Prediction Interval")+
  theme(plot.title = element_text(hjust=0)) 

plot(response_true,ylim=c(-0.4,0.4),type="l")
lines(pred_interval[,1],col="blue")
lines(pred_interval[,2],col="blue")
response_true
pred_interval
