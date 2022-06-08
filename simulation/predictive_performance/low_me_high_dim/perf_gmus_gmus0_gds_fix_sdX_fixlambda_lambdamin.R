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

########################
#library(hdme)
create_example_data <- function(n, p, s = 5, sdX = 1, sdU = 0.5, 
                                sdEpsilon = 0.1, family = "gaussian") {
  # Independent true covariates with mean zero and standard deviation sdX
  X <- matrix(rnorm(n * p, sd = sdX), nrow = n, ncol = p)
  # if Gaussian, response has standard deviation sdEpsilon, and zero intercept
  # if binomial, response is binomial with mean (1 + exp(-X %*% beta))^(-1)
  beta <- c(-2, -1, 0.5, 1, 2, rep(0, p - s))
  
  if(family == "gaussian") {
    # True coefficient vector has s non-zero elements and p-s zero elements
    y <- X %*% beta + rnorm(n, sd = sdEpsilon)  
  } else if (family == "binomial") {
    # Need an amplification in the binomial case
    beta <- beta * 3
    y <- rbinom(n, size = 1, prob = (1 + exp(-X %*% beta))**(-1))
  }
  
  # The measurements W have mean X and standard deviation sdU. 
  # We assume uncorrelated measurement errors
  W <- X + matrix(rnorm(n * p, sd = sdU), nrow = n, ncol = p)
  
  return(list(X = X, W = W, y = y, beta = beta, sigmaUU = diag(p) * sdU))  
}

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
#######################################
#######################################
#######################################
#######################################
#######################################
#######################################
#######################################
#######################################



mos=log10(c(1.5,3,6,20,40,70,100,1e5,1e10,1e20)) #here it is: max(me)/sdX

est_beta_gmus<-c() ##   (5-fold-cv * B times)* 10 mos
est_beta_gmus0<-c() ##   (5-fold-cv * B times)* 10 mos
est_beta_gmus00<-c() ##   (5-fold-cv * B times)* 10 mos
est_beta_gds<-c()

cvloglik_gmus<-c()  ##   1cvloglik(5-fold-cv) 60 * B times)* 10 mos
cvloglik_gmus0<-c() ##   1cvloglik(5-fold-cv) 60 * B times)* 10 mos
cvloglik_gmus00<-c() ##   1cvloglik(5-fold-cv) 60 * B times)* 10 mos
cvloglik_gds<-c() ##   1cvloglik(5-fold-cv)60 * 30 rounds * B times)* 10 mos


performance<-c()
mean_perform_AUC_gmus<-c()
mean_perform_AUC_gds<-c()
perf_prop_nonzero_coef_in_true_nonzero_coef_gmusPlusgmus0<-c()
perf_prop_nonzero_coef_in_true_nonzero_coef_gds<-c()
counts_nr=0

for(MOS in 1:length(mos)){
  
  counts_nr=counts_nr+1
  cat(counts_nr,"-th me/sdx start \n")
  n=60
  p=100
  s=10
  #beta <- 2*c(-10, -10, 10, 2, 1, rep(0, p - s))
  beta <- 1*c(-20,18,-16,14, -12,10,-8,6, -4,2, rep(0, p - s))
  me<-mos[counts_nr]*c(rep(c(0.33,0.66, 1.00),times=4)[1:10],rep(1,p-s))
  sdU=diag(x=me,nrow=p)
  sdX=1
  B=50
  counts=0
  
  vec_gmus0_AUC<-c()
  vec_gmus_AUC<-c()
  vec_gmus00_AUC<-c()
  vec_gds_AUC<-c()
  prop_nonzero_coef_in_true_nonzero_coef_gmus<-c()
  prop_nonzero_coef_in_true_nonzero_coef_gmus0<-c()
  prop_nonzero_coef_in_true_nonzero_coef_gmus00<-c()
  prop_nonzero_coef_in_true_nonzero_coef_gds<-c()
  
  for (bstar in 1:B){ #B=50 datasets
    counts=counts+1
    set.seed(20211215+10*bstar)
    ll=generate_logistic_data(n=n,p=p,s=s,beta=beta,sdU=sdU,sdX=sdX)
    ##
    
    set.seed(11+bstar)
    ind<-c(rep(1,time=n/5),rep(2,time=n/5),rep(3,time=n/5),rep(4,time=n/5),rep(5,time=n/5))
    sample_ind<-sample(ind,replace = F)
    training_W<-ll$W[!(sample_ind==1),]
    test_W<-ll$W[sample_ind==1,]
    training_Y<-ll$y[!(sample_ind==1)]
    test_Y<-ll$y[sample_ind==1]
    
    dim(training_W)
    dim(test_W)
    length(training_Y)
    length(test_Y)
    
    pred_cv<-c()
    temp<-c()
    LAMBDA<-numeric(5) # fix lambda for each cv prediction, 5 CV => 5 lambdas
    
    for( CV in 1:5 ){ # 5-fold cross validation
      training_W<-ll$W[!(sample_ind==CV),]
      test_W<-ll$W[sample_ind==CV,]
      training_Y<-ll$y[!(sample_ind==CV)]
      test_Y<-ll$y[sample_ind==CV]
      
      merror=ll$sigmaUU
      gmus_fit <- gmus(training_W, training_Y,merror=merror, family = "binomial")
      count=0
      LAMBDA[CV]<-gmus_fit$lambda
      
      non_zero_coef<-apply(gmus_fit$beta,2,FUN = function(x) sum(abs(x)>1e-10) )
      temp_df<-data.frame(y=non_zero_coef,x=gmus_fit$delta)
      
      C0<-as.matrix(data.frame(x0=rep(1,nrow(temp_df)),x1=gmus_fit$delta,x1sq=gmus_fit$delta^2))
      d0<-temp_df$y
      A0<-matrix(c(0,0,-1),nrow=1)
      b0<-0
      
      fit_poly<-pracma::lsqlincon(C0,d0,A0,b0)
      
      
      a=fit_poly[3] #y=ax^2+bx+c
      #print(a)
      b=fit_poly[2]
      c=fit_poly[1]
      
      if (a>0 & (b*b-4*a*(c-s)<0) ){
        delta_esti=-b/2/a
        #print(delta_esti)
        if(delta_esti>=0) gmus_fit2<-gmus(training_W, training_Y,merror=merror,delta=delta_esti,lambda = LAMBDA[CV] ,family = "binomial")
        if(delta_esti<0) gmus_fit2<-gmus(training_W, training_Y,merror=merror,delta=0,lambda = LAMBDA[CV] ,family = "binomial")
      }
      
      
      if (a>0 & (b*b-4*a*(c-s)>=0) ){
        delta_esti=(-b-sqrt(b*b-4*a*(c-s)))/2/a
        #print(delta_esti)
        if(delta_esti>=0) gmus_fit2<-gmus(training_W, training_Y,merror=merror,delta=delta_esti,lambda = LAMBDA[CV] ,family = "binomial")
        if(delta_esti<0) gmus_fit2<-gmus(training_W, training_Y,merror=merror,delta=0,lambda = LAMBDA[CV] ,family = "binomial")
      }
      
      if (a==0){
        if (b==0) {
          fit_inv<-lm(y ~ I(1/x),data = temp_df[2:nrow(temp_df),])
          a=fit_inv$coefficients[1] # 
          b=fit_inv$coefficients[2] # y=a+b/x
          delta_esti=b/(s-a)
          if(delta_esti>=0) gmus_fit2<-gmus(training_W, training_Y,merror=merror,delta=delta_esti,lambda = LAMBDA[CV] ,family = "binomial")
          if(delta_esti<0) gmus_fit2<-gmus(training_W, training_Y,merror=merror,delta=0,lambda = LAMBDA[CV] ,family = "binomial")
        }
        if (b!=0){
          delta_esti=(s-c)/b
          #print(delta_esti)
          if(delta_esti>=0) gmus_fit2<-gmus(training_W, training_Y,merror=merror,delta=delta_esti,lambda = LAMBDA[CV] ,family = "binomial")
          if(delta_esti<0) gmus_fit2<-gmus(training_W, training_Y,merror=merror,delta=0,lambda = LAMBDA[CV] ,family = "binomial")
        }
      }  
      
      temp<-cbind(1/(1+exp(-test_W %*% gmus_fit2$beta)),test_Y)
      pred_cv<-rbind(pred_cv,temp)
      
      prop_nonzero_coef_in_true_nonzero_coef_gmus<-c(prop_nonzero_coef_in_true_nonzero_coef_gmus,
                                                     mean(order(abs(beta), decreasing = TRUE)[1:s] %in% order(abs(gmus_fit2$beta)[abs(gmus_fit2$beta)!=0], decreasing = TRUE)[1:s])) 
      #proportion: #(nonzero recovered coef)/#(true nonzero coef) 
      
      est_beta_gmus<-rbind(est_beta_gmus,t(gmus_fit2$beta))
      
    }  
    cvloglik_gmus<-c(cvloglik_gmus,pred_cv[,2]*log(pred_cv[,1])+(1-pred_cv[,2])*log(1-pred_cv[,1])) #cross-validated likelihood
    
    cat("gmus pred_cv has",nrow(pred_cv),"rows at",mos[counts_nr],"\n")
    
    ROCit_obj<-rocit(score=pred_cv[,1],class=as.numeric(pred_cv[,2]>0),method="bin")
    #plot(ROCit_obj)
    
    vec_gmus_AUC<-c(vec_gmus_AUC,ROCit_obj$AUC)
    
    
    #GMUS0 (modified GMUS with equal Measurement Error)
    pred_cv<-c()
    temp<-c()
    
    
    for( CV in 1:5 ){ # 5-fold cross validation
      training_W<-ll$W[!(sample_ind==CV),]
      test_W<-ll$W[sample_ind==CV,]
      training_Y<-ll$y[!(sample_ind==CV)]
      test_Y<-ll$y[sample_ind==CV]
      
      merror_equal=diag(rep(sqrt(mean(diag(ll$sigmaUU)^2)),times=length(diag(ll$sigmaUU))))
      gmus_fit <- gmus(training_W, training_Y,merror=merror_equal,lambda = LAMBDA[CV], family = "binomial")
      count=0
      
      non_zero_coef<-apply(gmus_fit$beta,2,FUN = function(x) sum(abs(x)>1e-10) )
      temp_df<-data.frame(y=non_zero_coef,x=gmus_fit$delta)
      
      C0<-as.matrix(data.frame(x0=rep(1,nrow(temp_df)),x1=gmus_fit$delta,x1sq=gmus_fit$delta^2))
      d0<-temp_df$y
      A0<-matrix(c(0,0,-1),nrow=1)
      b0<-0
      
      fit_poly<-pracma::lsqlincon(C0,d0,A0,b0)
      
      
      a=fit_poly[3] #y=ax^2+bx+c
      #print(a)
      b=fit_poly[2]
      c=fit_poly[1]
      
      if (a>0 & (b*b-4*a*(c-s)<0) ){
        delta_esti=-b/2/a
        #print(delta_esti)
        if(delta_esti>=0) gmus_fit2<-gmus(training_W, training_Y,merror=merror_equal,delta=delta_esti,lambda = LAMBDA[CV] ,family = "binomial")
        if(delta_esti<0) gmus_fit2<-gmus(training_W, training_Y,merror=merror_equal,delta=0,lambda = LAMBDA[CV] ,family = "binomial")
      }
      
      
      if (a>0 & (b*b-4*a*(c-s)>=0) ){
        delta_esti=(-b-sqrt(b*b-4*a*(c-s)))/2/a
        #print(delta_esti)
        if(delta_esti>=0) gmus_fit2<-gmus(training_W, training_Y,merror=merror_equal,delta=delta_esti,lambda = LAMBDA[CV] ,family = "binomial")
        if(delta_esti<0) gmus_fit2<-gmus(training_W, training_Y,merror=merror_equal,delta=0,lambda = LAMBDA[CV] ,family = "binomial")
      }
      
      if (a==0){
        if (b==0) {
          fit_inv<-lm(y ~ I(1/x),data = temp_df[2:nrow(temp_df),])
          a=fit_inv$coefficients[1] # 
          b=fit_inv$coefficients[2] # y=a+b/x
          delta_esti=b/(s-a)
          if(delta_esti>=0) gmus_fit2<-gmus(training_W, training_Y,merror=merror_equal,delta=delta_esti,lambda = LAMBDA[CV] ,family = "binomial")
          if(delta_esti<0) gmus_fit2<-gmus(training_W, training_Y,merror=merror_equal,delta=0,lambda = LAMBDA[CV] ,family = "binomial")
        }
        if (b!=0){
          delta_esti=(s-c)/b
          #print(delta_esti)
          if(delta_esti>=0) gmus_fit2<-gmus(training_W, training_Y,merror=merror_equal,delta=delta_esti,lambda = LAMBDA[CV] ,family = "binomial")
          if(delta_esti<0) gmus_fit2<-gmus(training_W, training_Y,merror=merror_equal,delta=0,lambda = LAMBDA[CV] ,family = "binomial")
        }
      }  
      
      temp<-cbind(1/(1+exp(-test_W %*% gmus_fit2$beta)),test_Y)
      pred_cv<-rbind(pred_cv,temp)
      
      
      prop_nonzero_coef_in_true_nonzero_coef_gmus0<-c(prop_nonzero_coef_in_true_nonzero_coef_gmus0,
                                                      mean(order(abs(beta), decreasing = TRUE)[1:s] %in% order(abs(gmus_fit2$beta)[abs(gmus_fit2$beta)!=0], decreasing = TRUE)[1:s])) 
      #proportion: #(nonzero recovered coef)/#(true nonzero coef)
      
      est_beta_gmus0<-rbind(est_beta_gmus0,t(gmus_fit2$beta))
    }
    
    cvloglik_gmus0<-c(cvloglik_gmus0,pred_cv[,2]*log(pred_cv[,1])+(1-pred_cv[,2])*log(1-pred_cv[,1])) #cross-validated likelihood
    ROCit_obj<-rocit(score=pred_cv[,1],class=as.numeric(pred_cv[,2]>0),method="bin")
    #plot(ROCit_obj)
    
    cat("gmus0 pred_cv has",nrow(pred_cv),"rows at",mos[counts_nr],"\n")
    
    vec_gmus0_AUC<-c(vec_gmus0_AUC,ROCit_obj$AUC)
    
    
    ####################################################################################
    
    pred_cv<-c()
    temp<-c()
    #gmus00 original GMUS
    for( CV in 1:5 ){ # 5-fold cross validation
      training_W<-ll$W[!(sample_ind==CV),]
      test_W<-ll$W[sample_ind==CV,]
      training_Y<-ll$y[!(sample_ind==CV)]
      test_Y<-ll$y[sample_ind==CV]
      
  
      gmus_fit <- hdme::gmus(training_W, training_Y,lambda = LAMBDA[CV], family = "binomial")
      count=0
      
      non_zero_coef<-apply(gmus_fit$beta,2,FUN = function(x) sum(abs(x)>1e-10) )
      temp_df<-data.frame(y=non_zero_coef,x=gmus_fit$delta)
      
      C0<-as.matrix(data.frame(x0=rep(1,nrow(temp_df)),x1=gmus_fit$delta,x1sq=gmus_fit$delta^2))
      d0<-temp_df$y
      A0<-matrix(c(0,0,-1),nrow=1)
      b0<-0
      
      fit_poly<-pracma::lsqlincon(C0,d0,A0,b0)
      
      
      a=fit_poly[3] #y=ax^2+bx+c
      #print(a)
      b=fit_poly[2]
      c=fit_poly[1]
      
      if (a>0 & (b*b-4*a*(c-s)<0) ){
        delta_esti=-b/2/a
        #print(delta_esti)
        if(delta_esti>=0) gmus_fit2<-hdme::gmus(training_W, training_Y,delta=delta_esti,lambda = LAMBDA[CV] ,family = "binomial")
        if(delta_esti<0) gmus_fit2<-hdme::gmus(training_W, training_Y,delta=0,lambda = LAMBDA[CV] ,family = "binomial")
      }
      
      
      if (a>0 & (b*b-4*a*(c-s)>=0) ){
        delta_esti=(-b-sqrt(b*b-4*a*(c-s)))/2/a
        #print(delta_esti)
        if(delta_esti>=0) gmus_fit2<-hdme::gmus(training_W, training_Y,delta=delta_esti,lambda = LAMBDA[CV] ,family = "binomial")
        if(delta_esti<0) gmus_fit2<-hdme::gmus(training_W, training_Y,delta=0,lambda = LAMBDA[CV] ,family = "binomial")
      }
      
      if (a==0){
        if (b==0) {
          fit_inv<-lm(y ~ I(1/x),data = temp_df[2:nrow(temp_df),])
          a=fit_inv$coefficients[1] # 
          b=fit_inv$coefficients[2] # y=a+b/x
          delta_esti=b/(s-a)
          if(delta_esti>=0) gmus_fit2<-hdme::gmus(training_W, training_Y,delta=delta_esti,lambda = LAMBDA[CV] ,family = "binomial")
          if(delta_esti<0) gmus_fit2<-hdme::gmus(training_W, training_Y,delta=0,lambda = LAMBDA[CV] ,family = "binomial")
        }
        if (b!=0){
          delta_esti=(s-c)/b
          #print(delta_esti)
          if(delta_esti>=0) gmus_fit2<-hdme::gmus(training_W, training_Y,delta=delta_esti,lambda = LAMBDA[CV] ,family = "binomial")
          if(delta_esti<0) gmus_fit2<-hdme::gmus(training_W, training_Y,delta=0,lambda = LAMBDA[CV] ,family = "binomial")
        }
      }  
      
      temp<-cbind(1/(1+exp(-test_W %*% gmus_fit2$beta)),test_Y)
      pred_cv<-rbind(pred_cv,temp)
      
      
      prop_nonzero_coef_in_true_nonzero_coef_gmus00<-c(prop_nonzero_coef_in_true_nonzero_coef_gmus00,
                                                      mean(order(abs(beta), decreasing = TRUE)[1:s] %in% order(abs(gmus_fit2$beta)[abs(gmus_fit2$beta)!=0], decreasing = TRUE)[1:s])) 
      #proportion: #(nonzero recovered coef)/#(true nonzero coef)
      
      est_beta_gmus00<-rbind(est_beta_gmus00,t(gmus_fit2$beta))
    }
    
    cvloglik_gmus00<-c(cvloglik_gmus00,pred_cv[,2]*log(pred_cv[,1])+(1-pred_cv[,2])*log(1-pred_cv[,1])) #cross-validated likelihood
    ROCit_obj<-rocit(score=pred_cv[,1],class=as.numeric(pred_cv[,2]>0),method="bin")
    #plot(ROCit_obj)
    
    cat("gmus0 pred_cv has",nrow(pred_cv),"rows at",mos[counts_nr],"\n")
    
    vec_gmus00_AUC<-c(vec_gmus00_AUC,ROCit_obj$AUC)
    
    #GDS
    
    AUC_gds<-c()
    pred_cv<-c()
    temp<-c()
    AUC_vec<-c()
    
    
    for (i in 1:5){
      training_W<-ll$W[!(sample_ind==i),]
      test_W<-ll$W[sample_ind==i,]
      training_Y<-ll$y[!(sample_ind==i)]
      test_Y<-ll$y[sample_ind==i]
      
      dim(training_W)
      dim(test_W)
      length(training_Y)
      length(test_Y)
        
      gds_fit<-hdme::gds(X = training_W,training_Y,lambda = LAMBDA[i],family = "binomial")
      #plot(gds_fit)
      temp<-cbind(1/(1+exp(-test_W %*% gds_fit$beta)),test_Y)
      pred_cv<-rbind(pred_cv,temp)
      prop_nonzero_coef_in_true_nonzero_coef_gds<-c(prop_nonzero_coef_in_true_nonzero_coef_gds,
                                                      mean(order(abs(beta), decreasing = TRUE)[1:s] %in% order(abs(gds_fit$beta)[abs(gds_fit$beta)!=0], decreasing = TRUE)[1:s]))
        #cat(i,"-th done \n")
        
      est_beta_gds<-rbind(est_beta_gds,t(gds_fit$beta))
    }
    
    ROCit_obj<-rocit(score=pred_cv[,1],class=as.numeric(pred_cv[,2]>0),method="bin")
    #plot(ROCit_obj)
    cvloglik_gds<-c(cvloglik_gds,pred_cv[,2]*log(pred_cv[,1])+(1-pred_cv[,2])*log(1-pred_cv[,1])) #cross-validated likelihood
    vec_gds_AUC<-c(vec_gds_AUC,ROCit_obj$AUC)
    
    cat(counts,"-th simulation done \n")
  }
  
  
  mean_perform_AUC_gmus<-c(mean_perform_AUC_gmus,mean(vec_gmus_AUC,na.rm=T))
  mean_perform_AUC_gds<-c(mean_perform_AUC_gds,mean(vec_gds_AUC,na.rm=T))
  
  
  #mean(vec_gmus_AUC>vec_gds_AUC,na.rm=T)
  performance<-cbind(performance,cbind(vec_gmus_AUC,vec_gmus0_AUC,vec_gmus00_AUC,vec_gds_AUC))
  
  perf_prop_nonzero_coef_in_true_nonzero_coef_gmusPlusgmus0<-cbind(perf_prop_nonzero_coef_in_true_nonzero_coef_gmusPlusgmus0,
                                                                   cbind(prop_nonzero_coef_in_true_nonzero_coef_gmus,
                                                                         prop_nonzero_coef_in_true_nonzero_coef_gmus0,
                                                                         prop_nonzero_coef_in_true_nonzero_coef_gmus00))
  perf_prop_nonzero_coef_in_true_nonzero_coef_gds<-cbind(perf_prop_nonzero_coef_in_true_nonzero_coef_gds,
                                                         prop_nonzero_coef_in_true_nonzero_coef_gds)
  
}



plot(y=mean_perform_AUC_gmus-mean_perform_AUC_gds,x=mos,type="l")

dim(est_beta_gmus)
dim(est_beta_gmus0)
dim(est_beta_gds)

length(cvloglik_gmus)
length(cvloglik_gmus0)
length(cvloglik_gds)

###### save data
#AUC
write.csv(performance,file="AUC_Gmus_Gmus0_Gmus00_Gds_fixsdX.csv")

write.csv(perf_prop_nonzero_coef_in_true_nonzero_coef_gmusPlusgmus0,
          "perf_prop_nonzero_coef_in_true_nonzero_coef_gmusPlusgmus0_fixsdX.csv")
write.csv(perf_prop_nonzero_coef_in_true_nonzero_coef_gds,
          "perf_prop_nonzero_coef_in_true_nonzero_coef_gds_fixsdX")


write.csv(est_beta_gmus,"est_beta_gmus_fixsdX.csv")
write.csv(est_beta_gmus0,"est_beta_gmus0_fixsdX.csv")
write.csv(est_beta_gmus00,"est_beta_gmus00_fixsdX.csv")
write.csv(est_beta_gds,"est_beta_gds_fixsdX.csv")


write.csv(cvloglik_gmus,"cvloglik_gmus_fixsdX.csv")
write.csv(cvloglik_gmus0,"cvloglik_gmus0_fixsdX.csv")
write.csv(cvloglik_gmus00,"cvloglik_gmus00_fixsdX.csv")
write.csv(cvloglik_gds,"cvloglik_gds_fixsdX.csv")






