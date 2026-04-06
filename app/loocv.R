loocv <- function(rho,phi,temp,y,l,nn){
  err_squ=rep(0,l)
  n <- dim(temp)[2]/nn
  for (k in 1:l){
    if (k==1){
      x <- c()
      for (i in 1:nn)
        x <- cbind(x, t(temp[,((i-1)*n+1):(i*n)])%*%as.matrix(phi[1,,i])/res)
      }else{
        x <- c()
        for (i in 1:nn)
          x <- cbind(x, t(temp[,((i-1)*n+1):(i*n)])%*%t(phi[1:k,,i])/res)
        }
    H_hat=x%*%solve(t(x)%*%x)%*%t(x)
    err_squ[k]=mean(((diag(rep(1,n))-H_hat)%*%y/
      (rep(1,n)-diag(H_hat)))^2)
  }
  return(which(err_squ==min(err_squ)))
}
