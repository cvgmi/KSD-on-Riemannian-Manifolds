# data should be a 3-dimensional array, i.e., a group of samples
# weight should be a function

library("MASS") # for the generalized inverse function ginv() 
library("rstiefel") # for rmf.matrix()
library("matrixcalc") # for is.positive.semi.definite()


KSD <- function(X, model){

  # Set default values in model
  model <- prep(X,model)
  
  # Initialize the matrix G^w_p:=( k^w_p(x_i,x_j) )_{i j}
  Gwp <- matrix(0,model$n,model$n)
  
  # Compute the matrix G^w_p:=( k^w_p(x_i,x_j) )_{i j}
  n <- model$n
  weight <- model$weight
  for (i in 1:n) for(j in 1:n) {
     x <- X[,,i]
     y <- X[,,j]
     G1 <- Grad_logp(x,model) + Grad_logk(x,y,model)
     G2 <- Grad_logp(y,model) + Grad_logk(y,x,model)
     kp <- square_kp(x,y,G1,G2,model) + constant_kp(x,y,model)
     Gwp[i,j] <- kp*weight(x)*weight(y)
    }
  
  # Compute the empirical KSDs U^w_n and V^w_n
  Vwn=mean(Gwp)
  Uwn=mean(Gwp[col(Gwp)!=row(Gwp)]) 
  
  # Print U^w_n and V^w_n
  
  print( list(Vwn = Vwn, Uwn = Uwn)  )
  return( list(Gwp = Gwp, Vwn = Vwn, Uwn = Uwn) )
  
}

MKSDE <- function(X,model,GoF=FALSE, n_prime=10000){

  # Set default values in model
  model <- prep(X,model)
  
  # Record the values
  s <- model$s
  N1 <- model$N1
  N2 <- model$N2
  n <- model$n
  weight <- model$weight
  
  # Set the default values of estimators
  F_u <- matrix(0,N1,N2)
  A_u <- matrix(0,N1,N1) 
  F_v <- F_u
  A_v <- A_u

  # Initialize the Q and b
  Q_u<- matrix(0,s,s)
  b_u<- rep(0,s)
  Q_v <- Q_u
  b_v <- b_u
  Q_ite <- Q_v
  b_ite <- b_v
  

  # Compute Q and b
  for ( i in 1:n ) for( j in 1:n ){
    x <- X[,,i]
    y <- X[,,j]
    Q_ite <- Q(x,y,model)*weight(x)*weight(y)
    b_ite <- b(x,y,model)*weight(x)*weight(y)
    Q_v <- Q_v + Q_ite
    b_v <- b_v + b_ite
    if(i!=j) {
      Q_u <- Q_u + Q_ite
      b_u <- b_u + b_ite
      }
    }
  
  # Compute the estimator
  # ginv() is the generalized inverse in case Q is singular

  theta_v <- - ginv(Q_v) %*% b_v
  theta_u <- - ginv(Q_u) %*% b_u
  
  if(is.positive.semi.definite(sym(Q_u))==FALSE)
    print("Q_u for U-stat is not positive semi-definite")
  
  
  if (model$family == "MF") {
    F_v <- matrix(theta_v,N1,N2)
    F_u <- matrix(theta_u,N1,N2)
    
    result <- list( F_v = F_v, F_u = F_u)
    
    if(GoF == TRUE){
    model$F <- F_v
    invisible(capture.output(sol_v <- KSD(X,model)))
    model$F <- F_u
    invisible(capture.output(sol_u <- KSD(X,model)))
    p_v <- G_of_F(X,sol_v,W="V",n_prime)
    p_u <- G_of_F(X,sol_u,W="U",n_prime)
    
    result$Vwn <- sol_v$Vwn
    result$Uwn <- sol_u$Uwn
    result$p_u <- p_u
    result$p_v <- p_v
    }
  }
  else if(model$family=="MB"){
    A_v <- matrix(theta_v,N1,N1)
    A_u <- matrix(theta_u,N1,N1)
    
    result <- list( A_v = A_v, A_u = A_u)
    
    if(GoF == TRUE){
    model$A <- A_v
    invisible(capture.output(sol_v <- KSD(X,model)))
    model$A <- A_u
    invisible(capture.output(sol_u <- KSD(X,model)))
    p_v <- G_of_F(X,sol_v,W="V",n_prime)
    p_u <- G_of_F(X,sol_u,W="U",n_prime)
    
    result$Vwn <- sol_v$Vwn
    result$Uwn <- sol_u$Uwn
    result$p_u <- p_u
    result$p_v <- p_v
    
    } 
  }
  else if(model$family=="MFB"){
    P_v <- matrix(theta_v,N1,N1+N2)
    P_u <- matrix(theta_u,N1,N1+N2)
    
    A_v <- P_v(N1,1:N1)
    A_u <- P_u(N1,1:N1)
    F_v <- P_v(N1,(N1+1):(N1+N2))
    F_u <- P_u(N1,(N1+1):(N1+N2))
    
    result <- list( A_v = A_v, A_u = A_u, F_u = F_u, F_v = F_v )
    
    if(GoF == TRUE){
    model$A <- A_v
    model$F <- F_v
    invisible(capture.output(sol_v <- KSD(X,model)))
    model$A <- A_u
    model$F <- F_u
    invisible(capture.output(sol_u <- KSD(X,model)))
    p_v <- G_of_F(X,sol_v,W="V",n_prime)
    p_u <- G_of_F(X,sol_u,W="U",n_prime)
    
    result$Vwn <- sol_v$Vwn
    result$Uwn <- sol_u$Uwn
    result$p_u <- p_u
    result$p_v <- p_v
    
    }
  }
  
  return(result)
}


### Other Functions

# This function set default values in the model list
prep <- function(X,model){
  
  # Extract the dimensions of data
  dim <- datadim(X)
  model$N1 <- dim[1] # number of rows of the samples
  model$N2 <- dim[2] # number of cols of the samples
  model$n <- dim[3] # number of samples
  model$r <- dim[4]
  
  # Set default model to model the data
  if(is.null(model$space)) model$space <- "Stiefel"
  if(is.null(model$family)) model$family <- "MF"
  if(is.null(model$kernel)) model$kernel <- "Gaussian"
  
  # Set the default parameters of the density family
  if(is.null(model$F)) model$F <- matrix(0,model$N1,model$N2) 
  if(is.null(model$A)) model$A <- matrix(0,model$N1,model$N1) 
  
  # Set the default parameters in the kernel function
  if(is.null(model$tau)) model$tau <- 1 
  if(is.null(model$beta)) model$beta <- 1 
  if(is.null(model$gamma)) model$gamma <- 1 
  
  # Set default weight q(x)/w(x) to be 1
  if(is.null(model$weight)) model$weight <- function(x) 1
  
  # Set the dimension s of zeta
  if( model$family == "MF") model$s <- model$N1*model$N2 
  else if(model$family == "MB") model$s <- model$N1*model$N1 
  else if(model$family == "MFB") model$s <- model$N1*(model$N1+model$N2)
  
  return(model)
}

# Extract the dimension of the data
datadim <- function(x){
  # N1, num. of rows of the matrices
  # N2, num. of cols of the matrices
  # n, num. of samples(matrices)
  # r, rank of the matrices, which is the r in V_r(N) or G_r(N).
  dim <- dim(x)
  if(is.null(dim(x))) {
    # if data is a vector
    N1 <- 1 
    N2 <- 1 
    n <- length(x) 
    r <- 1
  } else if(length(dim)==2){
    N1 <- dim[1]
    N2 <- 1
    n <- dim[2]
    r <- 1
  }
  else if(length(dim)==3) {
    N1 <- dim[1] # number of rows of the samples
    N2 <- dim[2] # number of cols of the samples
    n <- dim[3] # number of samples
    r <- round(sum(x[,,1]*x[,,1]))
  }
  
  return(c(N1,N2,n,r))
  
}

# Generating Stiefel matrices from MF(F)
rMF <- function(n,F){
  X <- array( rep(F,n) , dim = c(nrow(F),ncol(F),n) )
  for(i in 1:n) X[,,i] <- rmf.matrix(F)
  return(X)
}

# Symmetrization and Skew-symmetrization of x
sym <- function(x) 0.5 * ( x + t(x) )
asym <- function(x) 0.5 * ( x - t(x) )

# Frobenius inner product of x,y
Frob <- function(x,y) sum( x * y )

# kernel functions including Gaussian and Inverse quadratic (InverseQ)
k <- function(x,y,model){
  if(model$kernel=="Gaussian") 
    exp( - model$tau/2 * norm(x-y,"F")^2 )
  else if( model$kernel == "InverseQ" ) 
    ( model$beta + norm(x-y,"F")^2 ) ^ { -model$gamma } 
}
  
# The Euclidean gradient of log p(x)
Grad_logp <- function(x,model) {
  # if (family=="MF"|family=="MB"|family=="MFB") 
   (model$A+t(model$A)) %*% x + model$F
}

# The Euclidean gradient of log k(x,y)
Grad_logk <- function(x,y,model){
  if( model$kernel == "Gaussian") model$tau*y
  else if( model$kernel == "InverseQ" ) 2*model$gamma*y/( model$beta + norm(x-y,"F")^2)
}

# The gradient of zeta(X)
Grad_zeta <- function(x,model){
  # X is a single matrix
  N1 <- model$N1
  N2 <- model$N2
  s <- model$s
  Grad <- array(0, dim = c(N1,N2,s))
  if( model$family == "MF" ) 
    Grad <- array(c(diag(N1*N2)), dim = c(N1,N2,s))
  else if( model$family == "MB" ){
    Gd <- array( c(diag(N1*N1)), dim = c(N1,N1,s))
    for ( i in 1:s) Grad[,,i] <- (Gd[,,i]+t(Gd[,,i])) %*% x
  }else if( model$family == "MFB"){
    Gd <- array(c(diag(N1*N1)),dim = c(N1,N1,N1*N1))
    for ( i in N1*N1 ) Grad[,,i] <- (Gd[,,i]+t(Gd[,,i])) %*% x
    Grad[ , , (N1*N1+1):s ] <- array(c(diag(N1*N2)),dim = c(N1,N2,N1*N2))
  }
  
  return(Grad)
}

# The gradient of eta(x)
Grad_eta <- function(x,model) 
# if (family=="MF"|family=="MB"|family=="MFB") 
{ 0 }

# which will also be used in the computation of MKSDE
square_kp <- function(x,y,G1,G2,model){
  # Here G means the gradient of everything
  if( model$space == "Stiefel" ){
    Frob( asym(G1 %*% t(x)), asym(G2 %*% t(y)) ) * k(x,y,model)
  }
  else if( model$space == "Grassmann" ){
   Frob( asym(sym(G1) %*% x ), asym(sym(G2) %*% y) ) * k(x,y,model)
  }
}

inside_kp <- function(x,G,model){
  # Here G means the gradient of everything
  if( model$space == "Stiefel" ) asym(G %*% t(x))
  else if( model$space == "Grassmann" ) asym(sym(G1)) %*% x 
}


constant_kp <- function(x,y,model){
  
  if(model$space=="Stiefel"){
    if (model$kernel=="Gaussian") 
      ConstantTerm <- model$tau/2 * ( model$N1 - 1 ) * Frob(x,y) * k(x,y,model)
    else if(model$kernel=="InverseQ"){
      ConstantTerm1 <- model$gamma * ( model$N1 - 1 ) * Frob(x,y) * ( model$beta + norm(x-y,"F")^2 )^{-model$gamma-1}
      ConstantTerm2 <- 4 * model$gamma * norm(asym(x%*%t(y)),"F")^2 * ( model$beta + norm(x-y,"F")^2)^{-model$gamma-2}
      ConstantTerm <- ConstantTerm1 + ConstantTerm2
    }
  }
  else if(model$space=="Grassmann"){
    if(model$kernel=="Gaussian") model$tau* ( model$N1 * Frob(x,y) - r^2 ) * k( x , y , model )
    else if(model$kernel=="InverseQ"){
      ConstantTerm1 <- 2 * model$gamma * ( model$N1 * Frob(x,y) - r^2 ) * ( model$beta + norm(x-y,"F")^2 )^{ - model$gamma - 1 }
      ConstantTerm2 <- 4 * model$gamma * norm(x-y,"F")^2 * ( model$beta + norm(x-y,"F")^2 )^{ - model$gamma - 2 }
      ConstantTerm <- ConstantTerm1 + ConstantTerm2
    }
  }
  
  return(ConstantTerm)
}

# The matrix A(X,Y)
Q <- function(x,y,model){

  # The gradient of zeta(X)
  
  G1 <- Grad_zeta(x,model)
  G2 <- Grad_zeta(y,model)
  
  Z1 <- matrix(0,model$N1*model$N1,model$s)
  Z2 <- Z1

  G1 <- Grad_zeta(x,model)
  G2 <- Grad_zeta(y,model)

  for(i in 1:model$s ){
    Z1[,i] <- c(inside_kp(x,G1[,,i],model))
    Z2[,i] <- c(inside_kp(y,G2[,,i],model))
  }

  Q <- t(Z2) %*% Z1 * k(x,y,model)

  
  # 
  # 
  # G1 <- Grad_zeta(x,model)
  # G2 <- Grad_zeta(y,model)
  # 
  # Q <- matrix(0,model$s,model$s)
  # for( i in 1:model$s ) for( j in 1:model$s )
  #    Q[i,j] <- square_kp(x,y,G1[,,i],G2[,,j],model)

  return(Q)
}

# The vector b(X,Y)
b <- function(x,y,model){
 
  # Grad_eta(x,model) + 
   G1 <- Grad_logk(x,y,model)
   G2 <- Grad_zeta(y,model)
   
   Z1 <- c(inside_kp(x,G1,model))
   Z2 <- matrix(0,model$N1*model$N1,model$s)

   for(i in 1:model$s ){
     Z2[,i] <- c(inside_kp(y,G2[,,i],model))
   }
   
   b <- t(Z2) %*% Z1 * k(x,y,model)
   
   # b <- rep(0,model$s)
   # for( i in 1:model$s) 
   #   b[i] <- square_kp(x,y,G1,G2[,,i],model)
   
   return(b)
}

G_of_F <- function(X, sol, W ="U", n_prime=10000){
  
  n <- nrow(sol$Gwp)
  
  # Compute the eigenvalues of the Gram matrix G^w_n:=
  lambda <- eigen(sol$Gwp/n)
  lambda <- lambda$values

  # Sample from the standard Gaussian
  Z <- rnorm(n*n_prime)
  Z <- matrix(Z,n,n_prime)
  
  # Compute the statistics and the distribution
  if(W=="U") {
    stat <- n*sol$Uwn
    gamma <- lambda %*% (Z^2-1)
  }
  else if(W=="V"){
    stat <- n*sol$Vwn
    gamma <- lambda %*% Z^2
  }
  
  # Compute the p_value
  p_value <- sum(stat<gamma)/n_prime
  return(p_value)
}


MLE <- function(X){
  dim <- datadim(X)
  n <- dim[3]
  r <- dim[4]
 return( rowSums(X,dims = 2)/n * r )
}

Comparison <- function(model,F){
  
  i <- 0
  
  sample_n <- seq(100,500,20)
  
  error <- matrix(0,3,length(sample_n))

  for( n in sample_n ) {
    
    X <- rMF(n, F)
    
    model <- prep(X,model)
    
    sol_MKSDE <- MKSDE(X,model)
    
    i <- i + 1
    
    F_hat_u <- sol_MKSDE$F_u
    F_hat_v <- sol_MKSDE$F_v
    F_MLE <- MLE(X)
    
    error[1,i] <- norm(F_hat_u-F,"F")
    error[2,i] <- norm(F_hat_v-F,"F")
    error[3,i] <- norm(F_MLE-F,"F")

  } 
  
  return(error)
    
}


plotErr <- function(error){
  
  sample_n <- seq(100,500,20)
  
  xlab <- "Number of samples"
  ylab <- "Frobenius distance"
  
  Max <- max(error)
  Min <- min(error)
  
  par(xpd=TRUE)
  
  plot(sample_n,error[3,], col = "green", type = "l", xlab = xlab, ylab = ylab, ylim = c(Min,Max))
  lines(sample_n,error[1,], col = "red" )
  lines(sample_n,error[2,], col = "blue" )
  
  legend(x = "top", inset = c(0,-0.2), bty = "n", horiz = TRUE, legend = c("MKSDE-U","MKSDE-V","MLE"), fill = c("red","blue","green"))
  
}


p_values <- function(model){
  
  F_0 <- matrix(c(1,1,1,0,0,0),3,2)
  n <- c(100,150,200,250,300)
  F <- list(0.3*F_0, F_0, 5*F_0)
  p_values_u <- matrix(0,length(n),length(F))
  p_values_v <- p_values_u
  

  for ( i in 1:length(n) ) for( j in 1:length(F) ){
      X <- rMF(n[i],F[[j]])
      sol <- MKSDE(X,model,GoF = TRUE)
      p_values_u[i,j] <- sol$p_u
      p_values_v[i,j] <- sol$p_v
  }
  
  return(list(p_values_u = p_values_u, p_values_v = p_values_v))
}


G <- matrix(0,1000,1000)