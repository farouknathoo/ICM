signal_sim <- function(P,K,T, n_E, n_M, beta_u=NULL, c_length, sigma2_a, sigma2_A, sigma2_u1,SNR){
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
  library(rgl)
  library(R.matlab)
  library(scatterplot3d)
  library(MASS)
  library(PottsUtils)
  library(MCMCpack)
  
  # Vertices data is obtained from the MATLAB as a mat file named "vert".
  path <- getwd()
  setwd(path)
  data <- readMat("./Data/vert.mat")
  
  vert.data <- data.frame(x = data$vert[,1],y = data$vert[,2],z = data$vert[,3])
  
  # Take the subset from all 5124 vertices.
  sub.index <- sample(1:5124,P,replace = F)
  sub.vert <- vert.data[sub.index,]
  r3dDefaults$windowRect <- c(0,45,780,667) 
  plot3d(sub.vert,pch=19,col="red")
  play3d( spin3d(rpm=3), duration=5)
  
  len.cubes <- c_length
  x.cut <- seq(-68,68+len.cubes,len.cubes)
  # Numeber of intervals on x axis
  n.x <- length(x.cut) - 1
  
  y.cut <- seq(-105,68+len.cubes,len.cubes)
  # Numeber of intervals on y axis
  n.y <- length(y.cut) -1
  
  z.cut <- seq(-49,78+len.cubes,len.cubes)
  # Numeber of intervals on z axis
  n.z <- length(z.cut) - 1

  # total number of voxels:
  n.v <- n.x*n.y*n.z
  
  # For each vertex, finding which intervals its x, y, z in. 
  vl <- cbind(findInterval(sub.vert$x,x.cut),findInterval(sub.vert$y,y.cut),findInterval(sub.vert$z,z.cut))

  # Mapping the indices into the labelling of each cube. 
  vert.Z <- rep(NA, P)
  for(i in 1:P)
  {
    vert.Z[i] <- vl[i,1] + (vl[i,2] -1)*n.x + (vl[i,3] -1)*(n.x*n.y)
  }
  
  # Connectivity Matrix 
  A <- diag(sample(seq(0.9,0.95,0.01),K-1 ))
  
  M <- matrix(0, nrow = n_M, ncol = T)
  E <- matrix(0, nrow = n_E, ncol = T)
  S <- matrix(0, nrow = P, ncol = T)
  
  # Simulate the time for mu(t) and mu_1(t) = 0 for all t. 
  mu <- matrix(0, nrow = K, ncol = T)
  
  par(mfrow=c(3,2)) 
  # when t = 1
  t <- 1 
  mu[2:K,1] <- mvrnorm(n = 1, mu = rep(0,K-1), diag(x = sigma2_u1*1, ncol = K-1, nrow = K-1))
  while (t < T)
  {
    mu[2:K,t+1] <- mvrnorm(n = 1, mu = A%*% mu[2:K,t], Sigma = diag(x = sigma2_a*1, ncol = K-1, nrow = K-1))
    t <- t+1
  }
  
  matplot(t(mu),type = "l",ylab="mu(t)",xlab = "t",lwd = 2)
  title(main = "Latent State Dynamics")
  
  if(is.null(beta_u)){
  # critical value for beta is 0.6313259 when K = 4
  beta_u <- 2/3*log(0.5*(sqrt(2) + sqrt(4*K - 2)))
  }
  # Set the variances of mixture components
  alpha <- c(apply(mu, 1, var))*0.10
  alpha[1] <- mean(alpha[2:4])
  
  # Forward Operator
  X_M  <- matrix(replicate(n_M*P,rnorm(1,0,1)), ncol=P,byrow = T)
  X_E  <- matrix(replicate(n_E*P,rnorm(1,0,1)), ncol=P,byrow = T)
  
  
  # Generate random samples from a Potts model
  # 3D array specifying cubes.
  mask <- array(1, dim = c(n.x, n.y, n.z))
  
  # Define the neighborhood structure(First order) and get the neighbor matrix.
  neiStruc <- matrix(c(2,2,0,0,
                       0,2,0,0,
                       0,0,0,0), nrow=3, byrow=TRUE)
  neighbors <- getNeighbors(mask, neiStruc)
  blocks <- getBlocks(mask, nblock=2)
  Z.state <- BlocksGibbs(1, nvertex = n.v,ncolor = K, neighbors = neighbors, blocks = blocks, beta = beta_u )
  
  
  # Map the vertices into cube and get the cube state for correspoding vertex.
  cube.state <- Z.state[vert.Z]
  r3dDefaults$windowRect <- c(0,45,780,667) 
  plot3d(sub.vert,col = cube.state)
  play3d( spin3d(rpm=3), duration=5)
  # Covariance Matrix
  H_M <- diag(1, nrow = n_M, ncol = n_M)
  H_E <- diag(1, nrow = n_E, ncol = n_E)
  
  # Without the noise
  for (t in 1:T)
  {
    #     # Measurement Noise
    #     eps_m  <- mvrnorm( n = 1, mu = rep(0,n_M), Sigma = diag(x = sigma2_M, ncol = n_M, nrow = n_M) )
    #     eps_e <- mvrnorm( n = 1, mu = rep(0,n_E), Sigma = diag(x = sigma2_E, ncol = n_M, nrow = n_M) ) 
    
    # Mapping the state to corresponding components.
    for (j in 1:1000)
    {
      # Mixture Components
      mix.comp <- c(rnorm(1, mu[1,t],sqrt(alpha[1])),
                    rnorm(1, mu[2,t],sqrt(alpha[2])),
                    rnorm(1, mu[3,t],sqrt(alpha[3])),
                    rnorm(1, mu[4,t],sqrt(alpha[4])))
      
      S[j,t] <- mix.comp[cube.state[j]]
    }  
    M[,t] <- X_M %*% S[,t] 
    E[,t] <- X_E %*% S[,t]
    
  }
  
  matplot(t(S),type = "l",col = cube.state,ylab="S(t)", xlab="t")
  title(main = "Time course of S(t)")
  
  matplot(t(M),type = "l",ylab="M(t)", xlab="t")
  title(main = "MEG signals without noise")
  
  matplot(t(E),type = "l",ylab="E(t)", xlab="t")
  title(main = "EEG signals without noise")
  
  
  # Add the measurement noise into the signals
  sigma2_M <- SNR*mean(apply(M, 2, var))
  sigma2_E <- SNR*mean(apply(E, 2, var))
  
  eps_m  <- mvrnorm( n = T, mu = rep(0,n_M), Sigma = diag(x = sigma2_M, ncol = n_M, nrow = n_M) )
  eps_e <- mvrnorm( n = T, mu = rep(0,n_E), Sigma = diag(x = sigma2_E, ncol = n_M, nrow = n_M) ) 
  
  M_noise <- M + eps_m
  E_noise <- E + eps_e
  
  
  matplot(t(M_noise),type = "l",ylab="M_nosie(t)", xlab="t")
  title(main = "MEG signals with noise")
  
  matplot(t(E_noise),type = "l",ylab="E_noise(t)", xlab="t")
  title(main = "EEG signals with noise")
  
  save(P,K,T,S,M_noise,alpha,E_noise,cube.state,X_M,X_E,H_E,H_M,beta_u,Z.state,blocks,neighbors,vert.Z,mu,file = "./Data/simulation.RData")  
}
 