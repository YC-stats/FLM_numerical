select_m=function(n,y,fpc_eigen,fpc_score){
  range_m=1:10 #upper bound can be changed
  mm=length(range_m)
  R_m=rep(0,mm)
  for (i in 1:mm){
    m=range_m[i]
    if (m==1){
      c_hat=mean(fpc_score[,1]*y)
      R_m[i]=-sum((c_hat/fpc_eigen[1:m])^2)+2/n/(n-1)*
        sum((fpc_score[,1]*y-c_hat)^2)/fpc_eigen[1]^2
    }else{
      c_hat=apply(sweep(fpc_score[,1:m],1,y,"*"),2,mean) #m*1 vector
      R_m[i]=-sum((c_hat/fpc_eigen[1:m])^2)+2/n/(n-1)*
        sum(apply((sweep(sweep(fpc_score[,1:m],1,y,"*"),
                       2,c_hat,"-"))^2,2,sum)/fpc_eigen[1:m]^2)
    }
  }
  m_star=#which.min(R_m)+1
    max(which.min(R_m),2)
}