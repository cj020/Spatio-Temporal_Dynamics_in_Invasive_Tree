library(raster)
library(fields)
library(RColorBrewer)
library(rasterVis)
library(mvtnorm)
library(gridExtra)
library(fBasics)
library(coda)
library(tidyverse)
library(egg)
library(Rcpp)
library(RcppArmadillo)
library(Matrix)

sourceCpp("/functions/PDE_fun.cpp")
sourceCpp("/functions/H_matrix.cpp")

# Zero-inflated binomial likelihood function.
myZIbinomlikelihood <- function(x, gamma, p, s, log = TRUE){
  output <- NA
  
  if(sum(x<0|x>s)>0){
    likelihd <- 0
  } else{
    likelihd <- I(x==0)*(gamma + (1-gamma)*(1-p)^s) +
      I(x>0)*((1-gamma)*dbinom(x, size = s, prob = p, log = FALSE))
  }
  
  if(log==TRUE){
    loglikelihd <- log(likelihd)
    output <- loglikelihd
  } else{
    output <- likelihd
  }
  
  return(output)
}

# The apply version.
loglikhd_apply_fun <- function(edge, resp, prob, gamma.v, size){
  llikhd <- ifelse(edge==0, 
                   myZIbinomlikelihood(x = resp * 100, 
                                       gamma = gamma.v,
                                       p = prob, 
                                       s = size,
                                       log = TRUE), 0)
  return(llikhd)
}

# Zero-inflated binomial sampler.
myrZIbinom <- function(n, gamma, p, s){
  output <- rep(0, n)
  
  for (i in 1:n) {
    r <- runif(1, 0, 1)
    
    if(r <= gamma){
      output[i] <- 0
    } else{
      output[i] <- rbinom(n = 1, size = s, prob = p)
    }
  }
  
  return(output)
}

# Neighborhood function: extract the neighborhood information, output a n by 4 matrix, where n=number of grids.
neighborhood <- function(raster){
  nn <- matrix(rep(0, length(raster[ ])), nrow = length(raster[]), ncol = 4)
  
  for(i in 1:dim(nn)[1]){
    loc <- adjacent(raster,i)[ ,2]
    
    ln <- loc[which((loc-dim(raster)[1])==i)]
    rn <- loc[which((loc+dim(raster)[1])==i)]
    bn <- loc[which((loc+1)==i)]
    tn <- loc[which((loc-1)==i)]
    
    nn[i,1] <- if(length(ln)>0){ln} else{0}
    nn[i,2] <- if(length(rn)>0){rn} else{0}
    nn[i,3] <- if(length(bn)>0){bn} else{0}
    nn[i,4] <- if(length(tn)>0){tn} else{0}
  }
  return(nn)
}

# 1st order propagator function: output the H matrix with dimension n by n, where n is the number of grids. 
propagator.1st <- function(neighbor, diffu.rate, 
                           growth.rate, growth.rate2, dx, dy, dt){
  H <- matrix(0, nrow = dim(neighbor)[1], ncol = dim(neighbor)[1])
  
  for (i in 1:dim(neighbor)[1]) {
    
    if(sum(neighbor[i,]>0)==4)
    {
      H[i,i] <- 1-2*diffu.rate[i]*dt/dx^2-2*diffu.rate[i]*dt/dy^2+
        growth.rate[i]*dt
      H[i,neighbor[i,1]] <- diffu.rate[i]*dt/dx^2
      H[i,neighbor[i,2]] <- diffu.rate[i]*dt/dx^2
      H[i,neighbor[i,3]] <- diffu.rate[i]*dt/dy^2
      H[i,neighbor[i,4]] <- diffu.rate[i]*dt/dy^2
    }
  }
  
  return(H)
}

# 2nd order propagator function for logistic growth, matrix version.
propagator.2nd <- function(neighbor, diffu.rate, 
                           growth.rate, growth.rate2, dx, dy, dt){
  H2 <- matrix(0, nrow = dim(neighbor)[1], ncol = dim(neighbor)[1])
  
  for (i in 1:dim(neighbor)[1]) {
    
    if(sum(neighbor[i,]>0)==4){
      H2[i,i] <- -growth.rate2[i]*dt
    }
  }
  
  return(H2)
}

# 2nd order propagator function for logistic growth, simplying to the vector version, since the matrix is diagonal.
propagator.2nd.vec <- function(neighbor, diffu.rate, 
                               growth.rate, growth.rate2, dx, dy, dt){
  H2.vec <- rep(0, dim(neighbor)[1])
  
  for (i in 1:dim(neighbor)[1]) {
    
    if(sum(neighbor[i,]>0)==4){
      H2.vec[i] <- -growth.rate2[i]*dt
    }
  }
  
  return(H2.vec)
}

# 2nd order propagator function for exponential growth, output a zero matrix.
propagator.2nd.exp.growth <- function(neighbor, diffu.rate, 
                                      growth.rate, growth.rate2, dx, dy, dt){
  H2 <- matrix(0, nrow = dim(neighbor)[1], ncol = dim(neighbor)[1])
  
  for (i in 1:dim(neighbor)[1]) {
    
    if(sum(neighbor[i,]>0)==4){
      H2[i,i] <- 0
    }
  }
  
  return(H2)
}

# 1st order propagator matrix for fokker planck equation.
propagator.FP.1st <- function(neighbor, diffu.rate, 
                              growth.rate, dx, dy, dt){
  H <- matrix(0, nrow = dim(neighbor)[1], ncol = dim(neighbor)[1])
  
  for (i in 1:dim(neighbor)[1]) {
    
    if(sum(neighbor[i,]>0)==4)
    {
      H[i,i] <- 1-2*diffu.rate[i]*dt/dx^2-2*diffu.rate[i]*dt/dy^2+
        growth.rate[i]*dt
      H[i,neighbor[i,1]] <- diffu.rate[neighbor[i,1]]*dt/dx^2
      H[i,neighbor[i,2]] <- diffu.rate[neighbor[i,2]]*dt/dx^2
      H[i,neighbor[i,3]] <- diffu.rate[neighbor[i,3]]*dt/dy^2
      H[i,neighbor[i,4]] <- diffu.rate[neighbor[i,4]]*dt/dy^2
    }
  }
  
  return(H)
}

# 2nd order propagator matrix for fokker planck equation with logistic growth.
propagator.FP.2nd <- function(neighbor, diffu.rate, 
                              growth.rate, dx, dy, dt){
  H2 <- matrix(0, nrow = dim(neighbor)[1], ncol = dim(neighbor)[1])
  
  for (i in 1:dim(neighbor)[1]) {
    
    if(sum(neighbor[i,]>0)==4){
      H2[i,i] <- -growth.rate[i]*dt
    }
  }
  
  return(H2)
}

# 2nd order propagator matrix for fokker planck equation with exponential growth, output a zero matrix.
propagator.FP.2nd.exp.growth <- function(neighbor, diffu.rate, 
                                         growth.rate, dx, dy, dt){
  H2 <- matrix(0, nrow = dim(neighbor)[1], ncol = dim(neighbor)[1])
  
  for (i in 1:dim(neighbor)[1]) {
    
    if(sum(neighbor[i,]>0)==4){
      H2[i,i] <- 0
    }
  }
  
  return(H2)
}

# Harmonic mean function: removing grid with 0, output=[mean(1/x)]^{-1}.
myharmonicmeanfun <- function(x, na.rm){
  result <- NA
  n <- length(x)
  
  if(sum(x==0) == n){
    result <- 0
  } else{
    vec <- x[x!=0]
    result <- 1/mean(1/vec)
  }
  
  return(result)
}

# Mean function: removing grid with +-Inf, output=mean(x).
mymeanfun <- function(x, na.rm){
  result <- NA
  n <- length(x)
  
  if(sum(x==Inf) + sum(x==-Inf) == n){
    result <- 0
  } else{
    vec <- x[x!=Inf & x!=-Inf]
    result <- mean(vec)
  }
  
  return(result)
}

# Division function with raster: output with rst format
myrstdividefun <- function(x, y, na.rm){
  n1 <- dim(x)[1]
  n2 <- dim(x)[2]
  result <- x[]/y[]
  one.mat <- matrix(rep(1, n1*n2), nrow = n1, ncol = n2)
  
  result.rst <- raster(one.mat, xmn = 0, xmx = 1, ymn = 0, ymx = 1, 
                       crs = NA)
  result.rst[] <- result
  
  return(result.rst)
}

# Division function with raster: x/y^2, output with rst format
myrstdivide2fun <- function(x, y, na.rm){
  n1 <- dim(x)[1]
  n2 <- dim(x)[2]
  result <- x[]/(y[]^2)
  one.mat <- matrix(rep(1, n1*n2), nrow = n1, ncol = n2)
  
  result.rst <- raster(one.mat, xmn = 0, xmx = 1, ymn = 0, ymx = 1, 
                       crs = NA)
  result.rst[] <- result
  
  return(result.rst)
}

# Numerical solution of PDE with diffusion and logistic growth with homogenization: using Rcpp
# $\frac{\partial}{\partial t}\mu=\frac{\partial^2}{\partial \bs^2}(D\mu)+(G\mu)(1-\mu)$
mynumericalPDEhomofun.Rcpp_mat <- function(n, epsilon, initial.vec, stream.vec,
                                           diffu.vec, growth.vec, 
                                           length_of_grid, time, step, delta_t, 
                                           min.value = 0, max.value = 1){
  
  # diffusion coefficient and its homogenization, D_bar
  diffu.rst <- origin.rst
  diffu.rst[] <- diffu.vec
  
  diffu.bar <- aggregate(diffu.rst, fact = epsilon, na.rm = TRUE, 
                         fun = myharmonicmeanfun)
  
  # growth coefficient and its homogenization, into G_1_bar and G_2_bar
  growth.rst <- origin.rst
  growth.rst[] <- growth.vec
  
  growth.bar <- diffu.bar * aggregate(myrstdividefun(growth.rst, diffu.rst), 
                                      fact = epsilon, 
                                      na.rm = TRUE, 
                                      fun = mymeanfun)
  
  growth2.bar <- diffu.bar * aggregate(myrstdivide2fun(growth.rst, diffu.rst), 
                                       fact = epsilon, 
                                       na.rm = TRUE, 
                                       fun = mymeanfun)
  
  # initial values and rescaling, c0=mean(mu*diffusion coefficient)
  mu0 <- initial.vec
  mu.rst <- origin.rst
  mu.rst[] <- mu0
  
  c0 <- rescaling.rst
  c0 <- aggregate(mu.rst * diffu.rst, fact = epsilon, na.rm = TRUE, fun = mymeanfun)

  cmat <- as.matrix(c0[])
  
  # input coef for the next functions
  diffu.rc <- diffu.bar[]
  growth.rc <- growth.bar[]
  growth2.rc <- growth2.bar[]
  
  # obtain the 1st order propagator in matrix (Rcpp)
  H1.mat <- propagator_1st_Rcpp(neighbor = neighbormatrix, 
                                diffu_rate = diffu.rc, 
                                growth_rate = growth.rc, 
                                growth_rate2 = growth2.rc,
                                dx = epsilon*length_of_grid, 
                                dy = epsilon*length_of_grid, 
                                dt = delta_t)
  
  # obtain the 2nd order propagator in vector (Rcpp)
  H2.vec <- propagator_2nd_vec_Rcpp(neighbor = neighbormatrix, 
                                    diffu_rate = diffu.rc, 
                                    growth_rate = growth.rc, 
                                    growth_rate2 = growth2.rc,
                                    dx = epsilon*length_of_grid, 
                                    dy = epsilon*length_of_grid, 
                                    dt = delta_t)
  
  # set up the result matrix, the grid with stream = 0
  cmat <- cmat * (1-stream.vec)
  c.rst <- rescaling.rst
  c.rst[] <- as.vector(cmat)
  
  # convert from c0 to mu0
  mu.rst <- disaggregate(c.rst, epsilon)/1
  mu.rst[] <- mu.rst[]/diffu.rst[]
  result <- matrix(mu.rst[], ncol = 1)
  result.c <- matrix(c.rst[], ncol = 1)
  
  # count step, to identify the 1st step
  step_count <- 1
  
  for (t in 2:(time*step+1)) {
    
    # update c, the grid with stream = 0
    cmat <- H1.mat %*% cmat + H2.vec * cmat^2
    cmat <- cmat * (1-stream.vec)
    
    if(step_count==step){ # only return the result at t=Step
      
      # convert c to mu
      c.rst <- rescaling.rst
      c.rst[] <- as.vector(cmat)
      mu.rst <- disaggregate(c.rst, epsilon)/1
      mu.rst[] <- mu.rst[]/diffu.rst[]
      
      # update result
      result <- cbind(result, mu.rst[])
      result.c <- cbind(result.c, c.rst[])
      
      step_count <- 1
    } else {
      step_count <- step_count+1
    }
  }
  
  return(list(result.mu = result, result.c = result.c, 
              diffu.bar = diffu.bar[], 
              growth.bar = growth.bar[],
              growth2.bar =growth2.bar[],
              H1 = H1.mat,
              H2 = H2.vec))
  
}

# Numerical solution of PDE with diffusion and logistic growth without homogenization:
# $\frac{\partial}{\partial t}\mu=\frac{\partial^2}{\partial \bs^2}(D\mu)+(G\mu)(1-\mu)$
mynumericalPDEfun <- function(n, initial.vec, stream.vec,
                              diffu.vec, growth.vec, 
                              length_of_grid, time, step, delta_t, 
                              min.value = 0, max.value = 1){
  
  # diffusion rate
  diffu.rc <- diffu.vec
  
  # growth rate
  growth.rc <- growth.vec
  
  # obtain 1st order propagator matrix
  H1.mat <- propagator.FP.1st(neighbor = neighbormatrix, 
                              diffu.rate = diffu.rc, 
                              growth.rate = growth.rc,
                              dx = length_of_grid, dy = length_of_grid, dt = delta_t)
  
  # obtain 2nd order propagator matrix
  H2.mat <- propagator.FP.2nd(neighbor = neighbormatrix, 
                              diffu.rate = diffu.rc, 
                              growth.rate = growth.rc,
                              dx = length_of_grid, dy = length_of_grid, dt = delta_t)
  
  # initial value
  mu0 <- initial.vec
  mumat <- as.matrix(mu0)
  
  # set up result 
  result <- matrix(mu0, ncol = 1)
  
  # count step, to identify the 1st step
  step_count <- 1
  
  for (t in 2:(time*step+1)) {
    
    # update mu, the grid with stream = 0
    mumat <- H1.mat %*% mumat + H2.mat %*% mumat^2
    mumat <- mumat * (1-stream.vec)
    
    if(step_count==step){ # only return the result at t=Step
      
      # obtain mu and update result
      mu <- as.vector(mumat)
      result <- cbind(result, mu)
      
      step_count <- 1
    } else {
      step_count <- step_count+1
    }
  }
  
  return(list(result.mu = result,
              diffu = diffu.rc, 
              growth = growth.rc,
              H1 = H1.mat,
              H2 = H2.mat))
}

# Numerical solution of PDE with diffusion and exponential growth without homogenization:
# $\frac{\partial}{\partial t}\mu=\frac{\partial^2}{\partial \bs^2}(D\mu)+G\mu$
mynumericalPDEfun_expgr <- function(n, initial.vec, stream.vec,
                                    diffu.vec, growth.vec, 
                                    length_of_grid, time, step, delta_t, 
                                    min.value = 0, max.value = 1){
  
  # diffusion rate
  diffu.rc <- diffu.vec
  
  # growth rate
  growth.rc <- growth.vec
  
  # obtain 1st order propagator matrix
  H1.mat <- propagator.FP.1st(neighbor = neighbormatrix, 
                              diffu.rate = diffu.rc, 
                              growth.rate = growth.rc,
                              dx = length_of_grid, dy = length_of_grid, dt = delta_t)
  
  # initial value for PDE
  mu0 <- initial.vec
  mumat <- as.matrix(mu0)
  
  # set up result 
  result <- matrix(mu0, ncol = 1)
  
  # count step, to identify the 1st step
  step_count <- 1
  
  for (t in 2:(time*step+1)) {
    
    # update mu, the grid with stream = 0
    mumat <- H1.mat %*% mumat
    mumat <- mumat * (1-stream.vec)
    
    if(step_count==step){ # only return the result at t=Step
      
      # obtain mu and update the result
      mu <- as.vector(mumat)
      result <- cbind(result, mu)
      
      step_count <- 1
    } else {
      step_count <- step_count+1
    }
  }
  
  return(list(result.mu = result,
              diffu = diffu.rc, 
              growth = growth.rc,
              H1 = H1.mat))
}

# Numerical solution of PDE with diffusion and exponential growth with homogenization: using Rcpp
# $\frac{\partial}{\partial t}\mu=\frac{\partial^2}{\partial \bs^2}(D\mu)+G\mu$
mynumericalPDEhomofun_expgr.Rcpp_mat <- function(n, epsilon, initial.vec, stream.vec,
                                                 diffu.vec, growth.vec, 
                                                 length_of_grid, time, step, delta_t, 
                                                 min.value = 0, max.value = 1){
  
  # diffusion coefficient and its homogenization, D_bar
  diffu.rst <- origin.rst
  diffu.rst[] <- diffu.vec
  
  diffu.bar <- aggregate(diffu.rst, fact = epsilon, na.rm = TRUE, 
                         fun = myharmonicmeanfun)
  
  # growth coefficient and its homogenization, G_1_bar, G_2_bar (G_2_bar will not be used)
  growth.rst <- origin.rst
  growth.rst[] <- growth.vec
  
  growth.bar <- diffu.bar * aggregate(myrstdividefun(growth.rst, diffu.rst), 
                                      fact = epsilon, 
                                      na.rm = TRUE, 
                                      fun = mymeanfun)
  
  growth2.bar <- diffu.bar * aggregate(myrstdivide2fun(growth.rst, diffu.rst), 
                                       fact = epsilon, 
                                       na.rm = TRUE, 
                                       fun = mymeanfun)
  
  # initial values and rescaling, c0=mean(mu*diffusion coefficient)
  mu0 <- initial.vec
  mu.rst <- origin.rst
  mu.rst[] <- mu0
  
  c0 <- rescaling.rst
  c0 <- aggregate(mu.rst * diffu.rst, fact = epsilon, na.rm = TRUE, fun = mymeanfun)

  cmat <- as.matrix(c0[])
  
  # rate coefficients for the next functions
  diffu.rc <- diffu.bar[]
  growth.rc <- growth.bar[]
  growth2.rc <- growth2.bar[]
  
  # obtain 1st order propagator matrix
  H1.mat <- propagator_1st_Rcpp(neighbor = neighbormatrix, 
                                diffu_rate = diffu.rc, 
                                growth_rate = growth.rc, 
                                growth_rate2 = growth2.rc,
                                dx = epsilon*length_of_grid, 
                                dy = epsilon*length_of_grid, 
                                dt = delta_t)
  
  # set up result
  cmat <- cmat * (1-stream.vec)
  c.rst <- rescaling.rst
  c.rst[] <- as.vector(cmat)
  
  # convert c0 to mu0
  mu.rst <- disaggregate(c.rst, epsilon)/1
  mu.rst[] <- mu.rst[]/diffu.rst[]
  
  result <- matrix(mu.rst[], ncol = 1)
  result.c <- matrix(c.rst[], ncol = 1)
  
  # count step, to identify the 1st step
  step_count <- 1
  
  for (t in 2:(time*step+1)) {
    
    # update c, the grid with stream = 0
    cmat <- H1.mat %*% cmat
    cmat <- cmat * (1-stream.vec)
    
    if(step_count==step){ # only return the result at t=Step
      
      # convert c to mu
      c.rst <- rescaling.rst
      c.rst[] <- as.vector(cmat)
      mu.rst <- disaggregate(c.rst, epsilon)/1
      mu.rst[] <- mu.rst[]/diffu.rst[]
      
      # update result
      result <- cbind(result, mu.rst[])
      result.c <- cbind(result.c, c.rst[])
      
      step_count <- 1
    } else {
      step_count <- step_count+1
    }
  }

  return(list(result.mu = result, result.c = result.c, 
              diffu.bar = diffu.bar[], 
              growth.bar = growth.bar[],
              H1 = H1.mat))
}