rm(list=ls())
source("gcv_b.R")
source("gcv_lamb.R")
source("gepsilon.R")
source("gAR0.R")
n=800
repli_b=1000
repli_s=1000
alpha=c(0.05,0.1)
res=100
a=11
d1=0.2
d2=0.5

t=(1:res)/res
kk=length(alpha)
v = 10# floor(n^(2/5))

b <- rep(0,a)
b[1]=0.8
b[2]=0.5
b[3]=-0.3
for (i in 4:a)
  b[i]=exp(-i)

basis=matrix(1,a,res)
basis[2:a,]=sqrt(2)*cos(pi*(1:(a-1))%*%t(t))
beta_true=as.vector(t(b)%*%basis) #res dimension
#cov_x <- c()

## Generate simulated data

# for (k in 1:repli_s){
#   set.seed(1678+n*(k-1))
#   u=matrix(runif(n*a,-sqrt(3),sqrt(3)),n,a)
#   xx=sweep(u,2,(1:a)^(-1),"*")
#   x_mat=xx%*%basis #n*res
#   eps=rnorm(n,0,1)
#   mean_x <- apply(x_mat,2,mean)
#   center_x <- sweep(x_mat,1,mean_x,"-")
#   cov_x <- cbind(cov_x,t(center_x)%*%center_x/n)
# }

#write.csv(cov_x, "cov.csv", row.names = FALSE)

## bootstrap for confidence bands

eig_val <- as.matrix(read.csv("eigenvalues_800.csv",header = FALSE))
eig_func <- as.matrix(read.csv("eigenfunctions_800.csv", header = FALSE))

count <- rep(0,2)
repli_b <- 1000
repli_s <- 1000
a_1 <- 1-1/sqrt(2)
a_2 <- 1+sqrt(2)

width = matrix(0,2,repli_s)


for (k in 1:repli_s){
  set.seed(1686 + k) 
  uu = matrix(runif(n*a,-sqrt(3),sqrt(3)),n,a)
  xx = sweep(uu,2,(1:a)^(-1),"*")
  eps = #rnorm(n,0,1)
    gepsilon(n,d1)
  #decay=exp(-0.5*(0:(a-1)))
  #xx=gAR0(n,d2,a,decay)  
  y=as.vector(xx%*%b+eps)
  
  x_mat = xx%*%basis #n*res
  #mean_x <- apply(x_mat,2,mean)
  #center_x <- sweep(x_mat,1,mean_x,"-")
  #cov_x <- cbind(cov_x,t(center_x)%*%center_x/n)
  
  rho <- eig_val[,k] # v*1
  eigen_phi <- t(eig_func[(v*(k-1)+1):(v*k),]) # len*v
  
  score <- x_mat%*%eigen_phi/res # n*v
  D <- diag(rho) # v*v
  
  lamb_hat <- gcv_lamb(y,score,v,D)$lamb
  ## estimation
  theta_hat <- solve(t(score)%*%score/n+lamb_hat*D)%*%t(score)%*%y/n # v*1
  beta_hat <- eigen_phi%*%theta_hat # len*1
  
  
  beta_boots <- matrix(0,res,repli_b)
  ## bootstrap procedure
  for (i in 1:repli_b){
    boots_seq <- a_1 + (a_2 - a_1)*rbinom(n, size = 1, prob = 1/3) # n*1
    #lamb <- gcv_b(y,score,D,boots_seq)$lamb
    M <- diag(boots_seq)
    beta_boots[,i] <- eigen_phi%*%(solve(t(score)%*%M%*%score/n+lamb_hat*D)%*%t(score)%*%M%*%y/n)
  }
  
  q <- apply(abs(sweep(beta_boots,1,beta_hat,"-")), 2,max) # repli_b*1
  ss <- sort(q)
  crit <- c(ss[floor(0.9*repli_b)], ss[floor(0.95*repli_b)])
  
  for (i in 1:2){
    if ( (sum( (beta_hat+crit[i] > beta_true) ) == res ) &
         (sum( (beta_hat-crit[i] < beta_true) ) == res ) ){
      count[i] <- count[i]+1
    }
    width[i,k] = 2*crit[i]
  }
}

print(round(apply(width,1,mean),3))
print(round(count/repli_s,3))
