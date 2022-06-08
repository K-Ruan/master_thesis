library(hdme)
library(ggplot2)
library(ROCit)
library(pracma)
library(quadprog)
### helper function

library(Rglpk)
library(glmnet)

# Logistic functions
logit <- function(x) (1+exp(-x))^(-1)
dlogit <- function(x) exp(-x)*(1+exp(-x))^(-2)

# Poisson functions
pois <- function(x) exp(x)
dpois <- function(x) exp(x)


###
musalgorithm <- function(W, y, merror, lambda, delta, weights = NULL){
  # We assume the first column of W is constants, i.e., intercept
  n <- dim(W)[1]
  p <- dim(W)[2]
  obj <- c(rep(1,p),rep(0,p)) #!there is no mistake here
  mat <- matrix(0,nrow=4*p, ncol=2*p)
  vec_merror<-diag(merror) #diagonal elements of the diagonal matrix
  
  # Weight matrix
  if(is.null(weights)){
    D <- diag(n)
  } else {
    D <- diag(weights)
  }
  
  #variable= (u_1,...,u_p,b_1,...,b_p)
  #obj_coef=(rep(0,p),rep(1,p))
  
  # Inequality constraint, -u_j - beta_j <= 0
  mat[1:p,1:p] <- -diag(p)
  mat[1:p,(p+1):(2*p)] <- -diag(p)
  
  # Inequality constraint, -u_j + beta_j <= 0
  mat[(p+1):(2*p),1:p] <- -diag(p)
  mat[(p+1):(2*p),(p+1):(2*p)] <- diag(p)
  
  transient_mat<-matrix(rep(c(1,vec_merror),p-1),nrow=p-1,ncol = p, byrow = T)
  
  # First "score function" constraint
  mat[(2*p+1),1:p] <- matrix(0, nrow=1, ncol=p) # intercept?
  mat[(2*p+2):(3*p),1:p] <- (-delta)*transient_mat
  mat[(2*p+1):(3*p),(p+1):(2*p)] <- 1/n*(t(W) %*% D %*% W)
  
  # Second "score function" constraint
  mat[(3*p+1),1:p] <- matrix(0, nrow=1, ncol=p)
  mat[(3*p+2):(4*p),1:p] <- (-delta)*transient_mat
  mat[(3*p+1):(4*p),(p+1):(2*p)] <- -1/n*(t(W) %*% D %*% W)
  
  rhs <- rep(0,(4*p))
  rhs[(2*p+1)] <- 1/n*(t(W[,1]) %*% D %*% y)
  rhs[(2*p+2):(3*p)] <- lambda + 1/n*(t(W[,-1]) %*% D %*% y)
  rhs[(3*p+1)] <- -1/n*(t(W[,1]) %*% D %*% y)
  rhs[(3*p+2):(4*p)] <- lambda - 1/n*(t(W[,-1]) %*% D %*% y)
  dir <- rep("<=",4*p)
  bounds <- list(lower=list(ind=1:(2*p), val=rep(-Inf,2*p)),
                 upper=list(ind=1:(2*p), val=rep(Inf,2*p)))
  
  
  
  
  bhat <- Rglpk::Rglpk_solve_LP(obj = obj, mat = mat, dir = dir,
                                rhs = rhs, bounds = bounds)$solution
  
  
  return(bhat[(p+1):(2*p)])
  
  
}

#merror:= covariate specific measurement error
#merror= diag(\sigma_1,...,\sigma_p) p*p matrix
mus_glm <- function(W, y, merror, lambda, delta, family = c("binomial", "poisson"), weights = NULL){
  
  family <- match.arg(family)
  
  if(family == "binomial") {
    mu <- logit
    dmu <- dlogit
  } else if(family == "poisson") {
    mu <- pois
    dmu <- dpois
  }
  
  
  n <- dim(W)[1]
  p <- dim(W)[2]
  
  bOld <- stats::rnorm(p)/p
  bNew <- stats::rnorm(p)/p
  IRLSeps <- 1e-7
  maxit <- 100
  count <- 1
  Diff1 <- 1
  Diff2 <- 1
  
  while(Diff1 > IRLSeps & Diff2 > IRLSeps & count < maxit){
    bOlder <- bOld
    bOld <- bNew
    V <- dmu(W%*%bOld)
    z <- W%*%bOld + (y - mu(W%*%bOld))/dmu(W%*%bOld)
    Wtilde <- c(sqrt(V)) * W
    ztilde <- c(sqrt(V)) * c(z)
    bNew <- musalgorithm(Wtilde, ztilde, merror, lambda, delta * sqrt(sum((V)^2)) / sqrt(n), weights)
    
    count <- count+1
    Diff1 <- sum(abs(bNew - bOld))
    Diff2 <- sum(abs(bNew - bOlder))
  }
  if(count >= maxit) print("Did not converge")
  return(bNew)
}

#merror:= covariate specific measurement error
#merror= diag(\sigma_1,...,\sigma_p) p*p matrix

gmus <- function(W, y, merror, lambda = NULL, delta = NULL,
                 family = "gaussian", weights = NULL) {
  
  family <- match.arg(family, choices = c("gaussian", "binomial", "poisson"))
  
  if(!is.null(weights) & length(weights) != nrow(W)) stop("weights vector must be one value per case")
  
  if(is.null(lambda)) {
    lambda <- glmnet::cv.glmnet(W, y, family = family)$lambda.min
  } else {
    stopifnot(all(lambda >= 0))
  }
  if(is.null(delta)) {
    delta <- seq(from = 0, to = 0.5, by = 0.02)
  } else {
    stopifnot(all(delta >= 0))
  }
  
  n <- dim(W)[1]
  p <- dim(W)[2] + 1
  W <- scale(W,center = T,scale = T)
  scales <- attr(W, "scaled:scale") 
  W <- cbind(rep(1,n), W) # add a column consisting of "1"s, representing the intercepts
  
  if(family == "gaussian") {
    fit <- sapply(delta, function(delta) musalgorithm(W, y, lambda, delta, weights))
  } else if(family %in% c("binomial", "poisson")) {
    fit <- sapply(delta, function(delta) mus_glm(W, y, merror, lambda, delta, family, weights))
  }
  
  
  fit <- list(intercept = fit[1, ],
              beta = matrix(fit[2:p, ] / scales, nrow = p - 1),
              family = family,
              delta = delta,
              lambda = lambda,
              num_non_zero = colSums(abs(fit[2:p, , drop = FALSE]) > 1e-10)
  )
  
  class(fit) <- "gmus"
  return(fit)
}

#########################################################
# data generating mechnism

generate_logistic_data<-function(n,p,s,beta,sdU,sdX){ #sdU: p*p diagonal matrix 
  vec_merror=diag(sdU) #length of vec == p
  dataXtrue<-c()
  dataXobs<-c()
  for(i in 1:p){
    onecol_trueX<-rnorm(n,sd=sdX)
    onecol_obsX<-onecol_trueX+rnorm(n,sd=vec_merror[i])
    dataXtrue<-c(dataXtrue,onecol_trueX)
    dataXobs<-c(dataXobs,onecol_obsX)
  }
  X<-matrix(data = dataXtrue,nrow = n, ncol = p, byrow = F)
  W<-matrix(data = dataXobs, nrow = n, ncol = p, byrow = F)
  y=rbinom(n, size = 1, prob = (1 + exp(-X %*% beta))**(-1))
  
  return(list(X = X, W = W, y = y, beta = beta, sigmaUU = sdU)) 
}

#######################################

