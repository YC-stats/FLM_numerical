##########multipier bootstrap#############
pboots_p<-function(S,x,f,t,m,c,time,R,basis,gg,p){
  SS=cbind(S[,m],(S[,(m+1):n]-S[,1:(n-m)]))/sqrt(m)
  r=dim(SS)[1]
  distr=rep(0,time)
  n=dim(x)[1]
  SS=cbind(S[,m],(S[,(m+1):n]-S[,1:(n-m)]))/sqrt(m)
  Sigma_inv=solve(t(x)%*%x/n+R)
  uboots=matrix(0,nrow=r,ncol=time)
  g=matrix(0,res,p)
  g_hat=matrix(0,res,p)
  Q_hat=array(0,c(res,time,p))
  M_r=rep(0,time)
  c_length=sum(c)
  ####addition####
  for (i in 1:time){
    u=rnorm(n-m+1,0,1)
    uboots[,i]=apply(sweep(SS,2,u,"*"),1,sum)/sqrt(n-m+1)
    }
  if (p==1){
    A_hat=sweep(basis[1:c[1],],1,f[1:c[1]],"/")
    Q=t(A_hat)%*%Sigma_inv%*%uboots
  }else if (p==2){
    A_hat=matrix(0,c_length,p*res)
    A_hat[1:c[1],1:res]=sweep(basis[1:c[1],],1,f[1:c[1]],"/")
    A_hat[(c[1]+1):c_length,(res+1):(2*res)]=
      sweep(basis[1:c[2],],1,f[(c[1]+1):(c[1]+c[2])],"/")
    Q=t(A_hat)%*%Sigma_inv%*%uboots
  }else{
    A_hat=matrix(0,c_length,p*res)
    A_hat[1:c[1],1:res]=sweep(basis[1:c[1],],1,f[1:c[1]],"/")
    A_hat[(c[1]+1):(c[1]+c[2]),(res+1):(2*res)]=
      sweep(basis[1:c[2],],1,f[(c[1]+1):(c[1]+c[2])],"/")
    A_hat[(c[1]+c[2]+1):c_length,(2*res+1):(3*res)]=
      sweep(basis[1:c[3],],1,f[(c[1]+c[2]+1):(c[1]+c[2]+c[3])],"/")
    Q=t(A_hat)%*%Sigma_inv%*%uboots
  }
  if (gg==1){
    M_r=apply(abs(Q),2,max)
    return(M_r)#vector
  }else{#function g is the std
      g=as.vector(apply(Q,1,sd))#(p*res)*1
      g_hat=rbind()
      mean_g=rep(0,p)
      for (i in 1:p){
        mean_g[i]=mean(g[((i-1)*res+1):(i*res)])
        g_hat=c(g_hat,g[((i-1)*res+1):(i*res)]/mean_g[i])
        #(g_hat[,i]<=0.1)
         #   g_hat[j,i]=max(g[((i-1)*res+1):(i*res)]/mean_g[i]/100)
        #}
        }
      M_r=apply(abs(sweep(Q,1,g_hat,"/")),2,max)
      g_mat=matrix(0,res,p)
      for (i in 1:p)
        g_mat[,i]=g_hat[((i-1)*res+1):(i*res)]
      list(M_r=M_r,g_hat=g_mat)
      }
  }