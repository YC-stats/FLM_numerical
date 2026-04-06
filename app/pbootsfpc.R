##########multipier bootstrap#############
pbootsfpc<-function(S,x,f,t,m,a,time,R,eig){
  SS=cbind(S[,m],(S[,(m+1):n]-S[,1:(n-m)]))/sqrt(m)
  r=dim(SS)[1]
  distr_1=rep(0,time)
  distr_2=rep(0,time)
  Sigma_inv=solve(t(x)%*%x/n+R)
  uboots1=matrix(0,nrow=r,ncol=time)
  for (i in 1:time){
    u=rnorm(n-m+1)
    uboots1[,i]=apply(sweep(SS,2,u,"*"),1,sum)/
      sqrt(n-m+1)
  }
  basis1=sweep(eig[[1]][,1:a[1]],2,f[[1]],"/")
  basis2=sweep(eig[[2]][,1:a[2]],2,f[[2]],"/")
  #basis3=sweep(eig[[3]][,1:a[3]],2,f[[3]],"/")
  g1=apply(basis1%*%Sigma_inv[1:a[1],]%*%uboots1,1,sd)#res*1
  g2=apply(basis2%*%Sigma_inv[(a[1]+1):(a[1]+a[2]),]%*%
             uboots1,1,sd)
  # g3=apply(basis3%*%Sigma_inv[(a[1]+a[2]+1):
  #           (a[1]+a[2]+a[3]),]%*%uboots1,1,sd)
  g1_hat=g1/mean(g1)
  g2_hat=g2/mean(g2)
  #g3_hat=g3#/mean(g3)
  uboots=matrix(0,nrow=r,ncol=time)  
  for (i in 1:time){
    u=rnorm(n-m+1)
    uboots=apply(sweep(SS,2,u,"*"),1,sum)/sqrt(n-m+1)
    Q1=as.vector(basis1%*%Sigma_inv[1:a[1],]%*%
                   uboots/g1_hat) #res*1 vector
    Q2=as.vector(basis2%*%Sigma_inv[(a[1]+1):(a[1]+a[2]),]%*%
                   uboots/g2_hat)
    # Q3=as.vector(basis3%*%Sigma_inv[(a[1]+a[2]+1):
    #              (a[1]+a[2]+a[3]),]%*%uboots/g3_hat)
    distr_1[i]=max(abs(Q1)) #,Q3
    distr_2[i]=max(abs(Q2))
  }
  list(distr_1=distr_1,distr_2=distr_2,g1=g1_hat,g2=g2_hat,#g3=g3_hat,
       g1_or=g1,g2_or=g2)#,g3_or=g3)
}