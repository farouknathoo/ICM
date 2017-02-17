# Simulation Study with all vertices used.

library(rgl)
library(R.matlab)
library(scatterplot3d)
library(MASS)
library(PottsUtils)
library(MCMCpack)
rm(list = ls())

P=3000
K=3
T=200
n_E=250
n_M=300
#c_length=20
#sigma2_a=1
#sigma2_A=1
#sigma2_u1=1
SNR=0.1
#beta_u=NULL

#this seems unstable so set it externally so that it works on your environment
#set it so that the current directory is ICM
#setwd(dir = "../../ICM/")

# The main working directory is always under "./ICM"
# R code for Bayesian source Reconstruction for MEG and EEG Data.
# P is the sample size of vertices.
# K is number of states.
# sigma2_a is the variance componenet of MVN used in autoregressive model for mu(t)
# sigma2_u1 is the variance component of MVN used for mu(t) when t = 1. 
# c_length is the length for the cube side.
# T is the time points. 
# n_E is the number of EEG sensors.
# n_M is the number of MEG sensors. 
# beta_U is the phase transition value. By default, it will be computed from beta_u <- 2/3*log(0.5*(sqrt(2) + sqrt(4*K - 2)))
# SNR is the measurement noise ratio. 

# Vertices data is obtained from the MATLAB as a mat file named "vert".

data <- readMat("./Data/vert.mat")
vert.data <- data.frame(x = data$vert[,1],y = data$vert[,2],z = data$vert[,3])
# sub.vert <- vert.data
# Take the subset from all 5124 vertices.
 #sub.index <- sample(1:5124,P,replace = F)
 #save(sub.index,file = "./Data/sub_index_3k.RData")
 load("./Data/sub_index_3k.RData")
 sub.vert <- vert.data[sub.index,]
#  r3dDefaults$windowRect <- c( 0,45,742,703)
P <- dim(sub.vert)[1]
#  play3d( spin3d(rpm=3), duration=5)

Dist<-matrix(0,nrow=P,ncol=P)
for(i in 1:P)
{
  x<-sub.vert$x[i]
  y<-sub.vert$y[i]
  z<-sub.vert$z[i]
  Dist[i,] <- sqrt( (sub.vert$x - x)^2 + (sub.vert$y - y)^2 + (sub.vert$z - z)^2)   
}

index1 <- 1834
c1 <- which( Dist[index1,] < 20)

which(sub.vert$x < -20 & sub.vert$y< -50 & sub.vert$z < -10)
index2 <-120
c2 <- which( Dist[index2,] < 20)

Z.state <- rep(0,P)
Z.state[c1] <- 2
Z.state[c2] <- 3
Z.state[-c(c1,c2)] <- 1

plot3d(sub.vert,pch=19,col=Z.state)

plot3d(sub.vert[-c(c1,c2),],pch=19,col="red")
points3d(sub.vert[c1,],pch=19,col="black")  
points3d(sub.vert[c2,],pch=19,col="blue")  
play3d( spin3d(rpm=3), duration=5)

t <- 1:T
s1 <- 20*dnorm(t,mean = 85,sd = 15)
#s1 <- 0.8*sin(0.1*t + 8) 
plot(t,s1,type='l')

s2 <- 20*dnorm(t,mean = 140,sd = 15)
points(t,s2,type = 'l',col="blue")

S <- matrix(0, nrow = P,ncol = T)
S[c1,] <- sweep(S[c1,] ,2, (-s1))
S[c2,] <- sweep(S[c2,] ,2, (-s2))
matplot(t(S),type = 'l',col = Z.state)

# Forward Operator
X_M  <- matrix(replicate(n_M*P,rnorm(1,0,sqrt(.1))), ncol=P,byrow = T)
X_E  <- matrix(replicate(n_E*P,rnorm(1,0,sqrt(.1))), ncol=P,byrow = T)

save(X_E,file = "./Data/X_E.RData")
save(X_M,file = "./Data/X_M.RData")

#Scale X for both simulated and real data
X_E <-  X_E / sqrt((1 / n_E)* sum(diag(X_E %*% t(X_E))))
X_M <-  X_M / sqrt((1 / n_M)* sum(diag(X_M %*% t(X_M))))


# len.cubes <- c_length
# x.cut <- seq(-68,68+len.cubes,len.cubes)
# # Numeber of intervals on x axis
# n.x <- length(x.cut) - 1
# 
# y.cut <- seq(-105,68+len.cubes,len.cubes)
# # Numeber of intervals on y axis
# n.y <- length(y.cut) -1
# 
# z.cut <- seq(-49,78+len.cubes,len.cubes)
# # Numeber of intervals on z axis
# n.z <- length(z.cut) - 1
# 
# # total number of voxels:
# n.v <- n.x*n.y*n.z
# 
# # For each vertex, finding which intervals its x, y, z in. 
# vl <- cbind(findInterval(sub.vert$x,x.cut),findInterval(sub.vert$y,y.cut),findInterval(sub.vert$z,z.cut))
# 
# # Mapping the indices into the labelling of each cube. 
# vert.Z <- rep(NA, P)
# for(i in 1:P)
# {
#   vert.Z[i] <- vl[i,1] + (vl[i,2] -1)*n.x + (vl[i,3] -1)*(n.x*n.y)
# }
# 
# # Connectivity Matrix 
# #YIN: REMOVE THIS AND HAVE A AS AN ARGUMENT TO THE FUNCTION
# A <- diag(sample(seq(0.9,0.95,0.01),K-1 ))

M <- matrix(0, nrow = n_M, ncol = T)
E <- matrix(0, nrow = n_E, ncol = T)
# S <- matrix(0, nrow = P, ncol = T)

# # Simulate the time for mu(t) and mu_1(t) = 0 for all t. 
# mu <- matrix(0, nrow = K, ncol = T)
# 
# par(mfrow=c(3,2)) 
# # when t = 1
# t <- 1 
# mu[2:K,1] <- mvrnorm(n = 1, mu = rep(0,K-1), diag(x = sigma2_u1*1, ncol = K-1, nrow = K-1))
# while (t < T)
# {
#   mu[2:K,t+1] <- mvrnorm(n = 1, mu = A%*% mu[2:K,t], Sigma = diag(x = sigma2_a*1, ncol = K-1, nrow = K-1))
#   t <- t+1
# }
# 
# matplot(t(mu),type = "l",ylab="mu(t)",xlab = "t",lwd = 2)
# title(main = "Latent State Dynamics")
# 
# if(is.null(beta_u)){
#   # critical value for beta is 0.6313259 when K = 4
#   beta_u <- 2/3*log(0.5*(sqrt(2) + sqrt(4*K - 2)))
# }
# # Set the variances of mixture components
# alpha <- c(apply(mu, 1, var))*0.10
# alpha[1] <- mean(alpha[2:4])
# 
# # Forward Operator
# X_M  <- matrix(replicate(n_M*P,rnorm(1,0,sqrt(.1))), ncol=P,byrow = T)
# X_E  <- matrix(replicate(n_E*P,rnorm(1,0,sqrt(.1))), ncol=P,byrow = T)
# 
# 
# # Generate random samples from a Potts model
# # 3D array specifying cubes.
# mask <- array(1, dim = c(n.x, n.y, n.z))
# 
# # Define the neighborhood structure(First order) and get the neighbor matrix.
# neiStruc <- matrix(c(2,2,0,0,
#                      0,2,0,0,
#                      0,0,0,0), nrow=3, byrow=TRUE)
# neighbors <- getNeighbors(mask, neiStruc)
# blocks <- getBlocks(mask, nblock=2)
# Z.state <- BlocksGibbs(1, nvertex = n.v,ncolor = K, neighbors = neighbors, blocks = blocks, beta = beta_u )
# 
# 
# # Map the vertices into cube and get the cube state for correspoding vertex.
# cube.state <- Z.state[vert.Z]
# r3dDefaults$windowRect <- c(0,45,780,667) 
# plot3d(sub.vert,col = cube.state)
# play3d( spin3d(rpm=3), duration=5)
# # Covariance Matrix
H_M <- diag(1, nrow = n_M, ncol = n_M)
H_E <- diag(1, nrow = n_E, ncol = n_E)

# Without the noise
for (t in 1:T)
{
  #   #     # Measurement Noise
  #   #     eps_m  <- mvrnorm( n = 1, mu = rep(0,n_M), Sigma = diag(x = sigma2_M, ncol = n_M, nrow = n_M) )
  #   #     eps_e <- mvrnorm( n = 1, mu = rep(0,n_E), Sigma = diag(x = sigma2_E, ncol = n_M, nrow = n_M) ) 
  #   
  #   # Mapping the state to corresponding components.
  #   for (j in 1:P)
  #   {
  #     # Mixture Components
  #     #YIN: NEED TO CHANGE THIS FOR GENERAL K
  #     mix.comp <- c(rnorm(1, mu[1,t],sqrt(alpha[1])),
  #                   rnorm(1, mu[2,t],sqrt(alpha[2])),
  #                   rnorm(1, mu[3,t],sqrt(alpha[3])),
  #                   rnorm(1, mu[4,t],sqrt(alpha[4])))
  #     
  #     S[j,t] <- mix.comp[cube.state[j]]
  #   }  
  M[,t] <- X_M %*% S[,t] 
  E[,t] <- X_E %*% S[,t]
  
}

matplot(t(S),type = "l",ylab="S(t)", xlab="t")
title(main = "Time course of S(t)")

matplot(t(M),type = "l",ylab="M(t)", xlab="t")
title(main = "MEG signals without noise")

matplot(t(E),type = "l",ylab="E(t)", xlab="t")
title(main = "EEG signals without noise")


# Add the measurement noise into the signals
sigma2_M <- SNR*mean(apply(M, 1, var))
sigma2_E <- SNR*mean(apply(E, 1, var))

eps_m  <- t(mvrnorm( n = T, mu = rep(0,n_M), Sigma = diag(x = sigma2_M, ncol = n_M, nrow = n_M) ))
eps_e <- t(mvrnorm( n = T, mu = rep(0,n_E), Sigma = diag(x = sigma2_E, ncol = n_E, nrow = n_E) ) )

M_noise <- M + eps_m
E_noise <- E + eps_e


matplot(t(M_noise),type = "l",ylab="M_nosie(t)", xlab="t")
title(main = "MEG signals with noise")

matplot(t(E_noise),type = "l",ylab="E_noise(t)", xlab="t")
title(main = "EEG signals with noise")

Y_M <- M_noise
Y_E <- E_noise

save(Y_M,file = "./Data/Y_M.RData")  
save(Y_E,file = "./Data/Y_E.RData")  
save(H_M,file = "./Data/H_M.RData")  
save(H_E,file = "./Data/H_E.RData")  
# save(X_E,file = "./Data/X_E.RData")
# save(X_M,file = "./Data/X_M.RData")
# save(cube.state,file = "./Data/cube_state.RData")
# save(A,file = "./Data/A.RData")
save(sigma2_E,file = "./Data/sigma2_E.RData")
save(sigma2_M,file = "./Data/sigma2_M.RData")
# save(sigma2_a,file = "./Data/sigma2_a.RData")
# save(mu,file = "./Data/mu.RData")
# save(sub.vert,file = "./Data/sub_vert.RData")
save(S, file = "./Data/S.RData")




# P=1000;K=4;T=200;n_E=300;n_M=250;c_length=20;sigma2_a=1;sigma2_A=1;sigma2_u1=1;SNR=0.1;beta_u=NULL

# signal_sim(P = 1000, n_E = 300,K = 4,T = 200,n_M = 250,beta_u = NULL,c_length = 20,sigma2_a = 1,sigma2_A = 1,sigma2_u1 = 1,SNR = 0.1)



