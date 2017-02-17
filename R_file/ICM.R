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
  library(nnet)
  # Load the simulation data.
  load("./Data/Y_E.RData")
  load("./Data/Y_M.RData")
  load("./Data/H_M.RData")
  load("./Data/H_E.RData")
  load("./Data/X_E.RData")
  load("./Data/X_M.RData")
  #load("./Data/sub_vert.RData")
  load("./Data/S.RData")
  #load("./Data/cube_state.RData")
  #load("./Data/A.RData")
  load("./Data/sigma2_E.RData")
  load("./Data/sigma2_M.RData")
  #load("./Data/sigma2_a.RData")
  #load("./Data/mu.RData")
  #load("./Data/beta.RData")
  #load("./Data/alpha.RData")
  data <- readMat("./Data/vert.mat")
  vert.data <- data.frame(x = data$vert[,1],y = data$vert[,2],z = data$vert[,3])
  # Take the subset from all 5124 vertices.
  # sub.index <- sample(1:5124,P,replace = F)
  # save(sub.index,file = "./Data/sub_index.RData")
  load("./Data/sub_index.RData")
  sub.vert <- vert.data[sub.index,]
  
  # Built-in data 
  n_M  <- dim(Y_M)[1] #NUMBER OF MEG SENSORS
  n_E  <- dim(Y_E)[1] #NUMBER OF EEG SENSORS
  T <-  dim(Y_M)[2] #NUMBER OF TIME POINTS
  P <- dim(sub.vert)[1] #NUMBER OF BRAIN LOCATIONS
  K <- 3 #NUMBER OF MIXTURE COMPONENTS (MESO-STATES)
  true.S <- S ; S <- NULL #STORE THE TRUE VALUES
  #true.mu <- mu; mu <- NULL
  #true.A <- A; A <- NULL
  true.sigma2_E <- sigma2_E; sigma2_E <- NULL
  true.sigma2_M <- sigma2_M; sigma2_M <- NULL
  #true.sigma2_a <- sigma2_a; sigma2_a <- NULL
  #true.cube.state <- cube.state; cube.state <- NULL
  #true.beta <- beta_u; beta_u <- NULL
  #true.alpha <- alpha; alpha <- NULL
  
  # Initial Number of voxels - this input will be modified
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
  #The following two lines are used for REAL but not simulated data
  #For simulated data the model generates Y_M and Y_E that are already scaled
  #Y_M <- Y_M / sqrt((1 / n_M)* sum(diag(Y_M %*% t(Y_M))))
  #Y_E <- Y_E / sqrt((1 / n_E)* sum(diag(Y_E %*% t(Y_E))))
  
  #Scale X for both simulated and real data
  X_E <-  X_E / sqrt((1 / n_E)* sum(diag(X_E %*% t(X_E))))
  X_M <-  X_M / sqrt((1 / n_M)* sum(diag(X_M %*% t(X_M))))
  
  
  
  ####################################################################
  ####################################################################
  ## Hyper-parameters set-up
  a_E <- 0.1; b_E <- 0.1
  a_M <- 0.1; b_M <- 0.1
  a_a <- 0.1; b_a <- 0.1
  a_alpha <- 0.1; b_alpha <- 0.1
  sigma2_A <- 0.25
  sigma2_mu1 <- 1 # check sensitivity
  beta_u <- 0.5*2/3*log(0.5*(sqrt(2) + sqrt(4*K - 2)))
  
  ###################################################################
  ###################################################################
  ## Initializing values
  beta_u_eps<-0.1
  # beta_u <- beta_u + beta_u_eps
  beta <- runif(1,0,beta_u)
  #beta <- true.beta
  
  A <- diag(sample(seq(0.9,0.95,0.01),K-1))
  #A <- true.A
  
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
  cube.state <- as.vector(cube.state)
  #cube.state <- true.cube.state
  v.index.nz <- as.numeric(names(table(vert.Z)))
  cube.state[-v.index.nz] <- 1
  Z.state <- cube.state[vert.Z]
  
  # For the variance components 
  sigma2_a <- b_a/(a_a+1) #initialize at prior mode
  #sigma2_a <- true.sigma2_a
  
  sigma2_E <- median(apply(Y_E,1,var))*0.1 #allocate 10% of median data variance to noise
  #sigma2_E <- true.sigma2_E
  sigma2_M <- median( apply(Y_M,1,var))*0.1
  #sigma2_M <- true.sigma2_M
  
  alpha <- rep(b_alpha/(a_alpha+1),K)
  #alpha <- true.alpha
  
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
  
  #mu <- true.mu
  
  #####get initial values for sources use either:
  # #Moore-Penrose Inverse
  # Y.start<-rbind(Y_E,Y_M)
  # X.start<-rbind(X_E,X_M)
  # X.start.mpinv<-ginv(X.start)
  # S <- X.start.mpinv%*%Y.start
  
  # #simulate from the assumed prior
  S<-matrix(rnorm(P*T,mean = 0,sd = 1),nrow=P,ncol=T)
  # S<-matrix(0,nrow=P,ncol=T)
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
  
  cor(c(true.S),c(S))
  
  
  ## ICM updating algorithm and intial set-up.
  R  <- 11 #NUMBER OF ICM ITERATIONS TO RUN
  
  ##STORE VALUES
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
  
  #THIS WILL BE A LARGE DATA STRUCTURE AND SHOULD ONLY BE USED FOR TESTING
  # S_star <- array(0,dim = c(P, T, R))
  # S_star[, , 1] <- S
  
  Z.nv <- matrix(0,nrow = n.v, ncol = R)
  Z.nv[,1] <- cube.state
  
  beta_star <- rep(0,R)
  beta_star[1] <- beta
  
  inv_H_M <- solve(H_M)
  inv_H_E <- solve(H_E)
  
  W_1j_C1 <- rep(0,P)
  W_1j_C2 <- rep(0,P)
  
  HX.M <- matrix(0,nrow = n_M,ncol = P)
  HX.E <- matrix(0,nrow = n_E,ncol = P)
  
  W2j_C1 <- matrix(0, nrow = T, ncol = P )
  W2j_C2 <- matrix(0, nrow = T, ncol = P )
  
  
  for ( j in 1:P)
  {
    HX.M[,j] <- inv_H_M %*% X_M[,j]
    HX.E[,j] <- inv_H_E %*% X_E[,j]
    
    W_1j_C1[j] <- crossprod(X_M[,j],HX.M[,j])
    W_1j_C2[j] <- crossprod(X_E[,j],HX.E[,j])
    
    for(t in 1:T){
      W2j_C1[t,j] <- -2 * crossprod(Y_M[,t], HX.M[,j])
      W2j_C2[t,j] <- -2 * crossprod(Y_E[,t], HX.E[,j])
    }
  }
  
  # Define function of -H(beta) and -H'(beta)
  
  H_beta <- function(beta,cube.state,n.v,neighbors)
  {
    test <- matrix(cube.state[neighbors],ncol = 6,byrow = FALSE)
    term1 <- 2*beta*sum(rowSums(test == cube.state,na.rm = TRUE))
    nn.count<-apply(test,1,tabulate,nbins=K)
    term2 <-sum(log(colSums(exp(2*beta*nn.count))))
    # Return the negative of H_beta 
    result <- -1*(term1 - term2)
    result
  }
  
  HP_beta <- function(beta,cube.state,n.v,neighbors)
  {
    test <- matrix(cube.state[neighbors],ncol = 6,byrow = FALSE)
    term1 <- 2*sum(rowSums(test == cube.state,na.rm = TRUE))
    nn.count<-apply(test,1,tabulate,nbins=K)
    term2.1 <- 1 / colSums(exp(2*beta*nn.count))
    term2.2 <- 2*colSums(exp(2*beta*nn.count) * nn.count)
    term2 <- sum(term2.1 * term2.2)
    result <- -1*(term1 - term2)
    result
  }
  
  r <- 1
  cor(c(true.S),c(S))
  
  while (r < R) {
    time.iter<-proc.time()
    # Update the sigma2_M
    a_M_star <- a_M + T*n_M / 2
    temp <- Y_M - X_M %*% S
    b_M_star <- 1/2 * sum( diag(crossprod(temp,inv_H_M%*%temp))) + b_M
    sigma2_M_star[r+1] <- b_M_star / (a_M_star + 1)
    sigma2_M <- b_M_star / (a_M_star + 1)
    #
    # Update the sigma2_E
    a_E_star <- a_E + T*n_E /2
    temp <- Y_E - X_E %*% S
    b_E_star <- 1/2 * sum( diag(crossprod(temp,inv_H_E%*%temp))) + b_E
    sigma2_E_star[r+1] <- b_E_star / (a_E_star + 1)
    sigma2_E <- b_E_star / (a_E_star + 1)
    #
    # Update the  sigma2_a
    a_a_star <- a_a + (T -1) * (K - 1) / 2
    temp <- mu[2:K, 2:T]- A%*%mu[2:K,1:T-1]
    b_a_star <- 1/2 * sum( diag( crossprod(temp))) + b_a
    sigma2_a_star[r+1] <- b_a_star / (a_a_star + 1)
    sigma2_a <- b_a_star / (a_a_star + 1)
    
    # 
    # Update the vec(A)
    sKr_t <- 0
    vc <- 0
    for (t in 2:T)
    {
      Kr_t <- kronecker( t(mu[2:K,t - 1]), diag(1, K-1))
      sKr_t <-  sKr_t + crossprod(Kr_t)#t(Kr_t) %*% Kr_t
      vc <- vc  + crossprod(mu[2:K,t],Kr_t)#t(mu[2:K,t]) %*% Kr_t
    }
    
    C_1 <- (1 / sigma2_A) * diag(1,(K-1)^2) + (1/sigma2_a) * sKr_t
    V_1 <- t( (1/sigma2_a) * vc %*% solve(C_1))
    A <- matrix(V_1, nrow = K-1)
    vec_A_star[,r+1] <- as.vector(A)
    
    # Update for each alpha_l, for each alpha_l, it's a inverse-gamma distribution.
    #  t.temp<-proc.time()
    a_alpha_l_star <- rep(0,K)
    b_alpha_l_star <- rep(0,K)
    for ( l in 1:K)
    {
      a_alpha_l_star[l] <- a_alpha + 1/2 * T * sum(Z.state == l)
      if (length(which(Z.state == l))==1)
      {
        b_alpha_l_star[l] <- b_alpha + 1/2 * sum ( (sweep(t(S[which(Z.state == l),]), 2, mu[l,]))^2 )
      }
      else
      {
        b_alpha_l_star[l] <- b_alpha + 1/2 * sum ( (sweep(S[which(Z.state == l),], 2, mu[l,]))^2 )
      }
      alpha[l] <- b_alpha_l_star[l] / ( a_alpha_l_star[l] + 1)
    }
    alpha_star[,r+1] <- alpha
    #  t.temp<-proc.time()-t.temp
    # 
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
    mu[2:K,1] <- M_1
    
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
    # SD_j <- apply(D_j, 1:2,sum)
    M_3 <- t( ( t(apply(STD_jT,2,sum)) + 1/sigma2_a*t(mu[2:K,T-1]) %*% t(A)) %*% inv_B_3)
    mu[2:K,T] <- M_3
    
    mu_star[,,r+1] <- mu
    
    #Update the ( S_j(1),  S_j(2), . . . ,  S_j(T)) for j = 1, 2, 3, ..., P.
    
    #OLD VERSION WITH NO VECTORIZATION
    #   
    #   for (j in 1:P){
    #     j=1
    #     W_1j <- 1/sigma2_M * t(X_M[,j]) %*%inv_H_M %*% X_M[,j] + 1/sigma2_E*t(X_E[,j]) %*% inv_H_E %*% X_E[,j] + 1/alpha[Z.state[j]]
    #     W_2j <- rep(0,T)
    #     for (t in 1:T){
    #       W_2j[t] <- 1/sigma2_M*( -2*t(Y_M[,t]) %*% inv_H_M %*% X_M[,j] + 2*t(X_M[,-j] %*% S[,t][-j]) %*% inv_H_M %*%X_M[,j] ) + 1/sigma2_E*(-2*t(Y_E[,t]) %*% inv_H_E %*% X_E[,j] + 2*t( X_E[,-j] %*% S[,t][-j]) %*% inv_H_E %*%X_E[,j])-2*sum(mu[,t]/alpha)
    #     }
    #     Sigma_S_j <- (1/as.numeric(W_1j))*diag(1,nrow = T,ncol = T)
    #     mu_sj <- -1/2*Sigma_S_j%*%W_2j
    #     S[j,] <- mu_sj
    #     t.temp<-proc.time()-t.temp
    #   }
    #   S_star[,,r+1] <- S
    
    # new version with vectorization
    W_1j <- 1/sigma2_M * W_1j_C1 + 1/sigma2_E*W_1j_C2+ 1/alpha[Z.state]
    for (j in 1:P)
    {
      # j = 1
      # W_2j <- 1/sigma2_M*( W2j_C1[,j] + 2*crossprod(X_M[,-j] %*% S[-j,], HX.M[,j]) ) + 1/sigma2_E*( W2j_C2[,j] + 2*crossprod(X_E[,-j] %*% S[-j,], HX.E[,j])) -2*rowSums( crossprod(mu, diag(1/alpha,nrow = K,ncol = K)))
      W_2j <- 1/sigma2_M*( W2j_C1[,j] + 2*crossprod(S[-j,], t(X_M[,-j])%*% HX.M[,j])) + 1/sigma2_E*( W2j_C2[,j] + 2*crossprod(S[-j,],t(X_E[,-j])%*% HX.E[,j])) -2*rowSums( crossprod(mu, diag(1/alpha,nrow = K,ncol = K)))
      Sigma_S_j <- (1/as.numeric(W_1j[j]))*diag(1,nrow = T,ncol = T)
      mu_sj <- -1/2* Sigma_S_j %*% W_2j
      S[j,] <- mu_sj
    }
    #S_star[,,r+1] <- S
    
    # 
    # #THIS NEEDS TO BE CHECKED CAREFULLY BY YIN - WE WERE TIRED WHEN WE FINISHED IT
    # # Update the labelling of Z.
    P.Z <- matrix(0, nrow = n.v,ncol = K)
    n_r <- rep(0,n.v)
    n_r.index <- as.numeric(names(table(vert.Z)))
    n_r.values <- as.vector(table(vert.Z))
    n_r[n_r.index] <- n_r.values
    for (v in v.index.nz)
    #for (v in 1:n.v)
    {
      for (h in 1:K)
      {
        log.term1 <- (-T*n_r[v]/2)*log(alpha[h])
        
        if (n_r[v] == 0)
        {
          log.term2 <- 0
        }
        else
        {
          if(n_r[v] ==1)
          {
            S.term2 <- t(S[which(vert.Z == v),])
          }
          else
          {
            S.term2 <- S[which(vert.Z == v),]
          }
          log.term2 <- -sum ((sweep(S.term2, 2, mu[h,]))^2)/ (2*alpha[h])
        }
        term3.neighbors <-neighbors[v,neighbors[v,] != (n.v+1)]
        log.term3 <- 4*beta*sum(cube.state[term3.neighbors] ==h)
        
        term4.temp <- neighbors[term3.neighbors,]
        term4.temp <- replace(term4.temp,term4.temp == v,n.v +1)
        nn.state <- matrix(cube.state[term4.temp],ncol = 6,byrow = FALSE)
        nn.count<- apply(nn.state,1,tabulate,nbins=K)
        nn.count[h,] <-nn.count[h,]+1
        log.term4<-sum(log(colSums(exp(2*beta*nn.count))))
        P.Z[v,h] <-  log.term1 + log.term2+ log.term3 - log.term4
        #P.Z[v,h] <-  log.term1+ log.term3 - log.term4
      }
      P.Z[v,] <-  P.Z[v,] -  max(P.Z[v,]) #ADDED FOR NUMRECIAL STABILITY
      P.Z[v,] <- exp(P.Z[v,])
      P.Z[v,] <-  P.Z[v,] / sum( P.Z[v,])
      cube.state[v] <- which.is.max(P.Z[v,]) #UPDATE THE STATE FOR VOXEL
      Z.state <- cube.state[vert.Z]
    }
    Z.nv[,r+1] <- cube.state #STORE UPDATED STATES FROM THIS SWEEP
    #  
    # 
    #UPDATE FOR BETA NEEDS TO BE ADDED
    beta.opt <-nlminb(start = beta, objective = H_beta,gradient = HP_beta, cube.state = cube.state,n.v= n.v, neighbors = neighbors,lower = 0,upper = beta_u)
    if(beta.opt$convergence != 0){
      cat("Warning, nlminb didn't converge for beta, iteration = ", r,"\n")
    }
    beta <- beta.opt$par
    beta_star[r+1] <- beta.opt$par
    #END ICM ITERATION
    time.iter<-proc.time()-time.iter
    cor.iter <- cor(c(true.S),c(S))
    rmse <- sqrt(sum( ( true.S - S)^2) /(P*T))
    cat("Iteration = ",r,"Time = ",time.iter[3]," seconds", "corr.iter =",cor.iter,"rmse = ",rmse, "\n")
    r <- r+1
  } 