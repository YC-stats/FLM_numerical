##########multipier bootstrap#############
pboots1<-function(S,x,f_hat,t,m,c,time,R,basis,g){
  distr=rep(0,time)
  if (c==1){
  n=length(x)
  SS=c(S[,m],(S[,(m+1):n]-S[,1:(n-m)]))/sqrt(m)
  Sigma_inv=1/(t(x)%*%x/n+R)
  }else{
    n=dim(x)[1]
    SS=cbind(S[,m],(S[,(m+1):n]-S[,1:(n-m)]))/sqrt(m)
    Sigma_inv=solve(t(x)%*%x/n+R)
  }
  uboots=matrix(0,nrow=c,ncol=time)
  if (g==1){
    ####addition####
    for (i in 1:time){
      u=rnorm(n-m+1,0,1)
      uboots[,i]=apply(sweep(SS,2,u,"*"),1,sum)/
          sqrt(n-m+1)
      }
    Q=t(sweep(basis[1:c,],1,f_hat,'/'))%*%
        Sigma_inv%*%uboots
    M_r=apply(abs(Q),2,max)
    return(M_r)
  }else{#function g is the std
    for (i in 1:time){
      u=rnorm(n-m+1,0,1)
      uboots[,i]=apply(sweep(SS,2,u,"*"),1,sum)/
          sqrt(n-m+1)
      }
    Q=t(sweep(basis[1:c,],1,f_hat,'/'))%*%
        Sigma_inv%*%uboots #res*time
    gg=as.vector(apply(Q,1,sd))
    g_stan=mean(gg)
    g_hat=gg/g_stan
    for (j in 1:length(gg)){
      if (g_hat[j]<=0.1)
        g_hat[j]=max(gg/g_stan/100)
    }
    M_r=apply(abs(sweep(Q,1,g_hat,'/')),2,max)
    list(M_r=M_r,g_hat=g_hat)
  }
}