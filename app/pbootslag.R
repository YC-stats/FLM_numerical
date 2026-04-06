##########multipier bootstrap#############
#for legendre bases
pbootslag<-function(S,x,f,t,m,a,time,R,leg,z=NA){ #add g_theta
  SS=cbind(S[,m],(S[,(m+1):n]-S[,1:(n-m)]))/sqrt(m)
  r=dim(SS)[1]
  l=length(a)
  distr=matrix(0,time,l)
  Sigma_inv=solve(t(x)%*%x/n+R)
  uboots=matrix(0,nrow=r,ncol=time)
  
  basis <- array(0, c(a[1],res,l))
  for (i in 1:l)
    basis[,,i] = sweep(leg[[i]][1:a[i],],1,f[,i],"/")
  
  for (i in 1:time){
    u=rnorm(n-m+1)
    uboots[,i]=apply(sweep(SS,2,u,"*"),1,sum)/
      sqrt(n-m+1)
  }
  
  g=sapply(1:l,function(i){ # functional coefficients
    apply(t(basis[,,i])%*%Sigma_inv[((i-1)*a[i]+1):(i*a[i]),]%*%uboots,1,sd)
  },simplify=TRUE) #res*l
  
  g_hat=sweep(g,2,apply(g,2,mean),"/")
  Q=matrix(0,res,l)
  
  if (!any(is.na(z))){
    g_theta <- apply((Sigma_inv[(sum(a)+1):r,,drop=FALSE]%*%uboots),1,sd)
    q_theta=matrix(0,time,r-sum(a))
  }else{
    g_theta <- c()
    q_theta <- c()
  }
  
  for (i in 1:time){
    u=rnorm(n-m+1)
    uboots=apply(sweep(SS,2,u,"*"),1,sum)/sqrt(n-m+1)
  
    for (k in 1:l){
      Q[,k]=as.vector(t(basis[,,k])%*%Sigma_inv[((k-1)*a[k]+1):(k*a[k]),]%*%
                     uboots/g_hat[,k]) #res*1 vector
    }
    distr[i,]=apply(abs(Q),2,max) 
    
    if (!any(is.na(z))){
      q_theta[i,] <- abs(as.vector(Sigma_inv[(sum(a)+1):r,,drop=FALSE]%*%uboots)/g_theta)
    }
  }
  list(distr=cbind(distr,q_theta), g=list(g_hat=g_hat,g_theta=g_theta), g_or=g)
}