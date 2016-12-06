## Algorithm Implementaion for ICM.
## ICM needs to find the mode of conditional distribution.
## Algorithm and pseudo code could be found in the Model_Setup document.
## Simulaton data could be read from the Data folder directory. 
## The main working directory is always under "./ICM"
setwd(dir = "~/ICM/")
rm(list = ls())
library(rgl)
library(R.matlab)
library(scatterplot3d)
library(MASS)
library(PottsUtils)
library(MCMCpack)
library(bigRR)
# Load the simulation data.
load("./Data/Y_E.RData")
load("./Data/Y_M.RData")
load("./Data/H_M.RData")
load("./Data/H_E.RData")
load("./Data/X_E.RData")
load("./Data/X_M.RData")
load("./Data/sub_vert.RData")
load("./Data/S.RData")
ls() 

# Built-in data 
n_M  <- dim(Y_M)[1]
n_E  <- dim(Y_E)[1]
T <-  dim(Y_M)[2]
P <- dim(sub.vert)[1]
K <- 4
true.S <- S

# Pre-set Number of cubes 
n.v <- 400

d.length <- ((range(sub.vert$x)[2] - range(sub.vert$x)[1])*(range(sub.vert$y)[2] - range(sub.vert$y)[1])*(range(sub.vert$z)[2] - range(sub.vert$z)[1]) / n.v)^(1/3)
len.cubes <- d.length 

x.cut <- seq(floor(range(sub.vert$x)[1]),range(sub.vert$x)[2]+len.cubes,len.cubes)
# Numeber of intervals on x axis
n.x <- length(x.cut) - 1

y.cut <- seq(floor(range(sub.vert$y)[1]),range(sub.vert$y)[2]+len.cubes,len.cubes)
# Numeber of intervals on y axis
n.y <- length(y.cut) -1

z.cut <- seq(floor(range(sub.vert$z)[1]),range(sub.vert$z)[2]+len.cubes,len.cubes)
# Numeber of intervals on z axis
n.z <- length(z.cut) - 1

# Actual total number of voxels:
n.v <- n.x*n.y*n.z
n.v

# For each vertex, finding which intervals its x, y, z in. 
vl <- cbind(findInterval(sub.vert$x,x.cut),findInterval(sub.vert$y,y.cut),findInterval(sub.vert$z,z.cut))

# Mapping the indices into the labelling of each cube. 
vert.Z <- rep(NA, P)
for(i in 1:P)
{
  vert.Z[i] <- vl[i,1] + (vl[i,2] -1)*n.x + (vl[i,3] -1)*(n.x*n.y)
}



####################################################################
####################################################################
## Scale the measurement of MEG and EEG.
#YIN:crossprod can speed this up
M <- Y_M / sqrt((1 / n_M)* sum(diag(Y_M %*% t(Y_M))))
E <- Y_E / sqrt((1 / n_E)* sum(diag(Y_E %*% t(Y_E))))

X_E <-  X_E / sqrt((1 / n_E)* sum(diag(X_E %*% t(X_E))))
X_M <-  X_M / sqrt((1 / n_M)* sum(diag(X_M %*% t(X_M))))

x.sd <- apply(rbind(X_E,X_M),2,sd)
X <- scale(rbind(X_E,X_M))

X_E <- X[1:n_E,]
X_M <- X[(n_E+1):(n_E+n_M),]


####################################################################
####################################################################
## Hyper-parameters set-up
a_E <- 0.1; b_E <- 0.1
a_M <- 0.1; b_M <- 0.1
a_a <- 0.1; b_a <- 0.1
a_alpha <- 0.1; b_alpha <- 0.1
sigma2_A <- 0.25
sigma2_mu1 <- 1 # check sensitivity
beta_u <- 2/3*log(0.5*(sqrt(2) + sqrt(4*K - 2)))

###################################################################
###################################################################
## Initializing values
beta_u_eps<-0.05
beta <- runif(1,beta_u-beta_u_eps,beta_u)
A <- diag(sample(seq(0.9,0.95,0.01),K-1))
# 3D array specifying cubes.
mask <- array(1, dim = c(n.x, n.y, n.z))

# Define the neighborhood structure(First order) and get the neighbor matrix.
neiStruc <- matrix(c(2,2,0,0,
                     0,2,0,0,
                     0,0,0,0), nrow=3, byrow=TRUE)
# Neighborhood matrix
neighbors <- getNeighbors(mask, neiStruc)

# Blocks for Chequeboard Updating
blocks <- getBlocks(mask, nblock=2)

# Get intial state for each vertex from correspoding cube. 
cube.state <- BlocksGibbs(1, nvertex = n.v,ncolor = K, neighbors = neighbors, blocks = blocks, beta = beta)
Z.state <- cube.state[vert.Z]


# For the variance components 
sigma2_a <- b_a/(a_a+1) #initialize at prior mode
sigma2_E <- median(apply(Y_E,1,var))*0.1 #allocate 10% of median data variance to noise
sigma2_M <- median( apply(Y_M,1,var))*0.1
alpha <- rep(b_alpha/(a_alpha+1),K)

# Initial values for mu
mu <- matrix(0, nrow = K, ncol = T)
# when t = 1
t <- 1 
mu[2:K,1] <- mvrnorm(n = 1, mu = rep(0,K-1), diag(x = sigma2_mu1*1, ncol = K-1, nrow = K-1))
while (t < T)
{
  mu[2:K,t+1] <- mvrnorm(n = 1, mu = A%*% mu[2:K,t], Sigma = diag(x = sigma2_a*1, ncol = K-1, nrow = K-1))
  t <- t+1
}

#####get initial values using ridge regression
# Initial value for S
Y.start<-rbind(Y_E,Y_M)
X.start<-rbind(X_E,X_M)
X.int<-matrix(1,nrow(Y.start),1)
S <- matrix(0, nrow = P, ncol = T)
for (t in 1: T)
{
  temp<-bigRR(y=Y.start[,t],X=X.int,Z=X.start,lambda=NULL,only.estimates=FALSE)
  S[,t] <- diag(1/x.sd,nrow = P,ncol = P)%*%temp$u
}


# for (t in 1:T)
# {
#   for (j in 1:P)
#   {
#     # Mixture Components
#     mix.comp <- c(rnorm(1, mu[1,t],sqrt(alpha[1])),
#                   rnorm(1, mu[2,t],sqrt(alpha[2])),
#                   rnorm(1, mu[3,t],sqrt(alpha[3])),
#                   rnorm(1, mu[4,t],sqrt(alpha[4])))
#     S[j,t] <- mix.comp[Z.state[j]]
#   }  
# }

## ICM updating algorithm and intial set-up.
R  <- 100
sigma2_M_star <- rep(0,R)
sigma2_M_star[1]  <- sigma2_M

sigma2_E_star <- rep(0,R)
sigma2_E_star[1]  <- sigma2_E

sigma2_a_star <- rep(0,R)
sigma2_a_star <- sigma2_a

vec_A_star <- matrix(0,nrow = (K-1)^2, ncol = R)
vec_A_star[,1] <- as.vector(A)

alpha_star <- matrix(0, nrow = K, ncol = R)
alpha_star[,1] <- alpha

mu_star <- array(0,dim = c(K,T,R))
mu_star[, , 1] <- mu

S_star <- array(0,dim = c(P, T, R))
S_star[, , 1] <- S

Z.nv <- matrix(0,nrow = n.v, ncol = R)
Z.nv[,1] <- cube.state

beta_star <- rep(0,R)
beta_star[1] <- beta

inv_H_M <- solve(H_M)
inv_H_E <- solve(H_E)
r <- 1

while (r < R) {
  
  # Update the sigma2_M 
  a_M <- a_M + T*n_M / 2
  b_M <- 1/2 * sum( diag(t(Y_M - X_M %*% S) %*% inv_H_M %*% (Y_M - X_M %*% S))) + b_M
  sigma2_M_star[r+1] <- b_M / (a_M + 1)
  
  # Update the sigma2_E
  a_E <- a_E + T*n_E /2
  b_E <- 1/2 *sum( diag(t(Y_M - X_M %*% S) %*% inv_H_E %*% (Y_M - X_M %*% S))) + b_E
  sigma2_E <- b_E / (a_E + 1)
  sigma2_E_star[r+1] <- b_E / (a_E + 1)
  
  # Update the  sigma2_a
  a_a <- a_a + (T -1) * (K - 1) / 2
  b_a <- b_a + 1/2 * sum( diag( t( mu[2:K, 2:T]- A%*%mu[2:K,1:T-1]) %*% (mu[2:K, 2:T]  - A %*%mu[2:K, 1:T-1])))
  sigma2_a <- b_a / (a_a + 1)
  sigma2_a_star[r+1] <- b_a / (a_a + 1)
  
  # Update the vec(A)
  sKr_t <- 0
  vc <- 0
  for (t in 2:T)
  {
    Kr_t <- kronecker( t(mu[2:K,t - 1]), diag(1, K-1))
    sKr_t <-  sKr_t + t(Kr_t) %*% Kr_t
    vc <- vc  + t(mu[2:K,t]) %*% Kr_t
  }
  
  C_1 <- 1 / sigma2_A * diag(1,(K-1)^2) + 1/sigma2_a * sKr_t
  V_1 <- t( 1/sigma2_a * vc %*% solve(C_1))
  A <- matrix(V_1, nrow = K-1)
  vec_A_star[,r+1] <- as.vector(A)
  
  # Update for each alpha_l, for each alpha_l, it's a inverse-gamma distribution.
  a_alpha_l <- rep(0,4)
  b_alpha_l <- rep(0,4)
  for ( l in 1:K)
  {
    a_alpha_l[l] <- a_alpha + 1/2 * T * sum(Z.state == l)
    b_alpha_l[l] <- b_alpha + 1/2 * sum ( (sweep(S[which(Z.state == l),], 2, mu[l,]))^2 )
    alpha[l] <- b_alpha_l[l] / ( a_alpha_l[l] + 1)
  }
  alpha_star[,r+1] <- alpha
  
  # Update mu_l(t=1) for all l =1, ... K.
  D_j <- array(0, dim = c(K-1, K-1,P))
  STD_j <- matrix(0,P,K-1)
  for (j in 1:P)
  {
    if (Z.state[j] != 1)
    {
      D_j[, , j][Z.state[j]-1,Z.state[j]-1] <- 1 / alpha[Z.state[j]] 
    }
    STD_j[j,] <- t(rep(S[j,1],K-1)) %*% D_j[, , j]  
  }
  
  SD_j <- apply(D_j, 1:2,sum)
  B_1 <- SD_j + 1/sigma2_a * t(A)%*%A + 1 / sigma2_mu1 * diag(1,nrow = K-1, ncol = K-1)
  M_1 <- t( ( t(apply(STD_j,2,sum)) + 1/sigma2_a*t(mu[2:K,2]) %*% A) %*% solve(B_1))
  mu[,1] <- rbind(0, M_1)
  
  # Update mu_l(t) for all l =2 , ... K when 1 < t < T,
  B_2 <- SD_j + 1 / sigma2_a * (t(A)%*%A + diag(1,K-1,K-1))
  inv_B_2 <- solve(B_2)
  STD_jt <- matrix(0,P,K-1)
  time_interval <- seq(2,T-1)
  for ( t in time_interval)
  {
    for (j in 1:P)
    {
      STD_jt[j,] <- t(rep(S[j,t],K-1)) %*% D_j[, , j]
    } 
    SSTD <- t(apply(STD_jt,2,sum))
    M_2_1 <- 1/sigma2_a*t(mu[2:K,t+1])%*%A
    M_2_2 <- 1/sigma2_a*t(mu[2:K,t-1]) %*% t(A) 
    M_2 <- t((SSTD + M_2_1 + M_2_2) %*% inv_B_2)
    mu[,t] <- rbind(0,M_2)
  }
  
  #Update mu_l(T) for all l=2, ..., K, when t = T.
  B_3 <- SD_j +1 / sigma2_a*diag(1,K-1,K-1)
  inv_B_3 <- solve(B_3)
  STD_jT <- matrix(0,P,K-1)
  for (j in 1:P)
  {
    STD_jT[j,] <- t(rep(S[j,T],K-1)) %*% D_j[, , j]  
  }
  SD_j <- apply(D_j, 1:2,sum)
  M_3 <- t( ( t(apply(STD_jT,2,sum)) + 1/sigma2_a*t(mu[2:K,T-1]) %*% t(A)) %*% inv_B_3)
  mu[,T] <- rbind(0, M_3)
  
  mu_star[,,r+1] <- mu
  
  #Update the ( S_j(1),  S_j(2), . . . ,  S_j(T)) for j = 1, 2, 3, ..., P.
  
  for (j in 1:P){
    W_1j <- 1/sigma2_M * t(X_M[,j]) %*%inv_H_M %*% X_M[,j] + 1/sigma2_E*t(X_E[,j]) %*% inv_H_E %*% X_E[,j] + alpha[Z.state[j]]
    W_2j <- rep(0,T)
    for (t in 1:T){
      W_2j[t] <- 1/sigma2_M*( -2*t(Y_M[,t]) %*% inv_H_M %*% X_M[,j] + 2*t(X_M[,-j] %*% S[,t][-j]) %*% inv_H_M %*%X_M[,j] ) + 1/sigma2_E*(-2*t(Y_E[,t]) %*% inv_H_E %*% X_E[,j] + 2*t( X_E[,-j] %*% S[,t][-j]) %*% inv_H_E %*%X_E[,j])
    }
    Sigma_S_j <- solve(diag(c(W_1j),nrow = T,ncol = T))
    mu_sj <- -1/2*Sigma_S_j%*%W_2j
    S[j,] <- mu_sj
  }
  S_star[,,r+1] <- S
  
  # Update the labelling of Z. 
  
}



