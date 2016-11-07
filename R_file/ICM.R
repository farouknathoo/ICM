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
# Load the simulation data.
load("./Data/Y_E.RData")
load("./Data/Y_M.RData")
load("./Data/H_M.RData")
load("./Data/H_E.RData")
load("./Data/X_E.RData")
load("./Data/X_M.RData")
load("./Data/sub_vert.RData")
ls() 

# Built-in data 
n_M  <- dim(Y_M)[1]
n_E  <- dim(Y_E)[1]
T <-  dim(Y_M)[2]
P <- dim(sub.vert)[1]
K <- 4

# Pre-set Number of cubes 
n.v <- 220

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

M <- Y_M / sqrt(1 / n_M* sum(diag(Y_M %*% t(Y_M))))
E <- Y_E / sqrt(1 / n_E* sum(diag(Y_E %*% t(Y_E))))

X_E <-  X_E / sqrt(1 / n_E* sum(diag(X_E %*% t(X_E))))
X_M <-  X_M / sqrt(1 / n_M* sum(diag(X_M %*% t(X_M))))


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

beta <- runif(1,0,beta_u)
A <- diag(sample(seq(0.9,0.95,0.01),K-1))
# 3D array specifying cubes.
mask <- array(1, dim = c(n.x, n.y, n.z))

# Define the neighborhood structure(First order) and get the neighbor matrix.
neiStruc <- matrix(c(2,2,0,0,
                     0,2,0,0,
                     0,0,0,0), nrow=3, byrow=TRUE)
neighbors <- getNeighbors(mask, neiStruc)
blocks <- getBlocks(mask, nblock=2)
# Get intial state for each vertex from correspoding cube. 
beta_u  <- 2/3*log(0.5*(sqrt(2) + sqrt(4*K - 2)))
cube.state <- BlocksGibbs(1, nvertex = n.v,ncolor = K, neighbors = neighbors, blocks = blocks, beta = beta)
Z.state <- cube.state[vert.Z]


# For the variance components 
sigma2_a <- 1
sigma2_E <- 100
sigma2_M <- 100
alpha <- c(1,1,1.5,2)

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

# Initial value for S
S <- matrix(0, nrow = P, ncol = T)
for (t in 1:T)
{
  for (j in 1:P)
  {
    # Mixture Components
    mix.comp <- c(rnorm(1, mu[1,t],sqrt(alpha[1])),
                  rnorm(1, mu[2,t],sqrt(alpha[2])),
                  rnorm(1, mu[3,t],sqrt(alpha[3])),
                  rnorm(1, mu[4,t],sqrt(alpha[4])))
    S[j,t] <- mix.comp[Z.state[j]]
  }  
}


