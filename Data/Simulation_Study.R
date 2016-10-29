# R code for Bayesian source Reconstruction for MEG and EEG Data.
library(rgl)
library(R.matlab)
library(scatterplot3d)
library(MASS)
library(PottsUtils)
library(MCMCpack)
# Task 1. How to map the vertices into voxels?
# Let's first plot the all the vertices in 3D first.
# Vertices data is obtained from the MATLAB as a mat file named "vert".
path <- getwd()
setwd(path)
data <- readMat("vert.mat")

# Plot vertices in 3D dimention. There are couple ways to do this
# Using scatterplot3d() function from scatterplot3d package without roration.
# or using plot3d from rgl function with the rotation included.

vert.data <- data.frame(x = data$vert[,1],y = data$vert[,2],z = data$vert[,3])

# Take the subset from all 5124 vertices.
sub.index <- sample(1:5124,1000,replace = F)
sub.vert <- vert.data[sub.index,]
plot3d(sub.vert,pch=19,col="red")
# points3d(sub.vert,pch=19,col="blue")
#play3d(spin3d(axis = c(0,1,0),rpm = 8),duration = 20)
range(sub.vert$x)
range(sub.vert$y)
range(sub.vert$z)

# Get the volume 
volume <-  sum(abs(range(vert.data$x)[1]) +abs(range(vert.data$x)[2]) )* sum(abs(range(vert.data$y)[1]) +abs(range(vert.data$y)[2]) )*sum(abs(range(vert.data$z)[1]) +abs(range(vert.data$z)[2]) )

# Set the length of interval.
len.cubes <- 25

# 
x.cut <- seq(-68,68+len.cubes,len.cubes)
x.cut

# Numeber of intervals on x axis
n.x <- length(x.cut) - 1

y.cut <- seq(-105,68+len.cubes,len.cubes)
y.cut

# Numeber of intervals on y axis
n.y <- length(y.cut) -1

z.cut <- seq(-49,78+len.cubes,len.cubes)
z.cut

# Numeber of intervals on z axis
n.z <- length(z.cut) - 1
n.z

# total number of voxels:
n.v <- n.x*n.y*n.z

# For each vertex, finding which intervals its x, y, z in. 
vl <- cbind(findInterval(sub.vert$x,x.cut),findInterval(sub.vert$y,y.cut),findInterval(sub.vert$z,z.cut))
vl

# Mapping the indices into the labelling of each cube. 
vert.Z <- rep(NA, 1000)
for(i in 1:1000)
  {
 vert.Z[i] <- vl[i,1] + (vl[i,2] -1)*n.x + (vl[i,3] -1)*(n.x*n.y)
  }
vert.Z
table(vert.Z)

######################################################################
######################################################################
# Simulate one realization from the model. 
# Set the parameters for Hyper Prior.
T <- 50
K <- 4
P <- 1000
n_M <- 200
n_E <- 200

# Weak Variance Prior
a_alpha <- 0.01; b_alpha <- 0.01
a_a <- 0.01; b_a <- 0.01
a_M <- 0.02; b_M <- 0.01
a_E <- 0.01; b_E <- 0.02

sigma2_u1 <- 1
sigma2_A <- 1

beta <- 2/3*log(0.5*(sqrt(2) + sqrt(4*K - 2)))
beta
# critical value for beta is 0.6313259 when K = 4


M <- matrix(0, nrow = n_M, ncol = T)
E <- matrix(0, nrow = n_E, ncol = T)
S <- matrix(0, nrow = P, ncol = T)

mu <- matrix(0, nrow = K, ncol = T)
# when t = 1
mu[2:K,1] <- mvrnorm(n = 1, mu = rep(0,K-1), diag(x = 1, ncol = K-1, nrow = K-1))

# Connectivity Matrix 
A <- diag(sample(seq(0.5,0.8,0.1),K-1 ))

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
Z.state <- BlocksGibbs(100, nvertex = 252,ncolor = 4, neighbors = neighbors, blocks = blocks, beta = 0.9)

# Variance Components

sigma2_M <- rinvgamma(1, a_M, b_M)
sigma2_E <- rinvgamma(1, a_E, b_E)
sigma2_a <- rinvgamma(1, a_a, b_a)

 alpha <- replicate(K,rinvgamma(1,a_alpha, b_alpha))

alpha <- c(1,1,1,1)

sigma2_E <- 1
sigma2_M <- 1
# Covariance Matrix
H_M <- diag(1, nrow = n_M, ncol = n_M)
H_E <- diag(1, nrow = n_E, ncol = n_E)

for (t in 1:T)
  {
    # Measurement Noise
    eps_m  <- mvrnorm( n = 1, mu = rep(0,n_M), Sigma = diag(x = sigma2_M, ncol = n_M, nrow = n_M) )
    eps_e <- mvrnorm( n = 1, mu = rep(0,n_E), Sigma = diag(x = sigma2_E, ncol = n_M, nrow = n_M) ) 
     
    # Mixture Components
    mix.comp <- c(rnorm(1, mu[,1],sqrt(alpha[1])),
                  rnorm(1, mu[,2],sqrt(alpha[2])),
                  rnorm(1, mu[,3],sqrt(alpha[3])),
                  rnorm(1, mu[,4],sqrt(alpha[4])))
    
    # Map the verrices into cube and get the cube state. Then, map this state to corresponding components.
    cube.state <- Z.state[vert.Z]
    S[,t] <- mix.comp[cube.state]
   
    M[,t] <- X_M %*% S[,t] 
    E[,t] <- X_E %*% S[,t]
    
   if (t < T)
      {
        mu[2:K,t+1] <- mvrnorm(n = 1, mu = A%*% mu[2:K,t], Sigma = diag(x = 1, ncol = K-1, nrow = K-1))
      }
  }



