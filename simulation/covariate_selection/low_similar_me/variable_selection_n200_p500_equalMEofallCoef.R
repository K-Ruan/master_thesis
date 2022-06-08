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


smoothing<-function(x){
  output<-numeric(length(x))
  for (i in 1:(length(x)-1)){
    if (x[i]<= min(x[-(1:i)])) output[i]<-x[i]
    if (x[i] > min(x[-(1:i)])) {
      x[i:(i+which.min(x[-(1:i)]))]<-mean(x[i:(i+which.min(x[-(1:i)]))])
      output[i:(i+which.min(x[-(1:i)]))]<-mean(x[i:(i+which.min(x[-(1:i)]))])
    } 
  }
  output[length(x)]<-x[length(x)]
  return(output)
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
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



#mos=log10(c(1.5,3,6,20,40,70,100,1e5,1e10,1e20)) #here it is: max(me)/sdX

mos=0.2 #here it is: max(me)/sdX

est_beta_gmus<-c() ##   B=100
est_beta_gmus0<-c() ##   B=100
est_beta_gds<-c() ## B=100

performance<-c()
counts_nr=0

nonzero_coef_gmus<-c() #to calculate MonteCarlo mean of #non-zero coef in the elbow plot
nonzero_coef_gmus0<-c()
nonzero_coef_gds<-c()

for(MOS in 1:length(mos)){
  
  counts_nr=counts_nr+1
  cat(counts_nr,"-th me/sdx start \n")
  
  n=200
  p=500
  s=10
  beta <- c(rep(1,times=s),rep(0,times=p-s))
  
  set.seed(1)
  me<-mos*c(abs(rnorm(s,mean=1,sd=0.2)),abs(rnorm(p-s,mean=1,sd = 0.2)) )
  sdU=diag(x=me,nrow=p)
  sdX=1
  B=100
  
  for (bstar in 1:B){
    
    set.seed(20211215+10*bstar)
    ll=generate_logistic_data(n=n,p=p,s=s,beta=beta,sdU=sdU,sdX=sdX)
    
    ## gmus
    merror=ll$sigmaUU
    gmus_fit <- gmus(ll$W, ll$y,merror=merror, family = "binomial")
    
    grid_nonzero_coef<-apply(gmus_fit$beta,MARGIN = 2, FUN=function(x){
      sum(abs(x)>1e-10)
    })
    nonzero_coef_gmus<-rbind(nonzero_coef_gmus, t(grid_nonzero_coef)   )
    
    #gmus0
    merror_equal=diag(rep(sqrt(mean(diag(ll$sigmaUU)^2)),times=length(diag(ll$sigmaUU))))
    gmus0_fit <- gmus(ll$W, ll$y,merror=merror_equal,
                      lambda=gmus_fit$lambda, family = "binomial")
    
    grid_nonzero_coef<-apply(gmus0_fit$beta,MARGIN = 2, FUN=function(x){
      sum(abs(x)>1e-10)
    })
    nonzero_coef_gmus0<-rbind(nonzero_coef_gmus0, t(grid_nonzero_coef)   )
    
    #gds is not needed here
  }
}   


##LOAD dataset here
load(file="nonzero_coef_gmus_me02_sample200_equalME0coef.RData")
load(file="nonzero_coef_gmus0_me02_sample200_equalME0coef.RData")
load(file="nonzero_coef_gds_me02_sample200_equalME0coef.RData")

delta_grid<-seq(0,0.5,by=0.02)

MonteCarlo_NCmean_gmus<-colMeans(nonzero_coef_gmus)
slope_gmus<-(MonteCarlo_NCmean_gmus[1:25]-MonteCarlo_NCmean_gmus[2:26])/(delta_grid[1:25]-delta_grid[2:26])
smooth_slope_gmus<-smoothing(slope_gmus)

MonteCarlo_NCmean_gmus0<-colMeans(nonzero_coef_gmus0)
slope_gmus0<-(MonteCarlo_NCmean_gmus0[1:25]-MonteCarlo_NCmean_gmus0[2:26])/(delta_grid[1:25]-delta_grid[2:26])
smooth_slope_gmus0<-smoothing(slope_gmus0)

#save(nonzero_coef_gmus,file="nonzero_coef_gmus_me02_sample200_equalME0coef.RData")
#save(nonzero_coef_gmus0,file="nonzero_coef_gmus0_me02_sample200_equalME0coef.RData")
#save(nonzero_coef_gds,file="nonzero_coef_gds_me02_sample200_equalME0coef.RData")
########## first run till here


#then choose delta
plot(0:25,MonteCarlo_NCmean_gmus) # NC =  nonzero coef
plot(1:25,y=smooth_slope_gmus)

plot(0:25,MonteCarlo_NCmean_gmus0)
plot(1:25,y=smooth_slope_gmus0)


####### plot

par(mfrow=c(1,1))
plot(delta_grid,MonteCarlo_NCmean_gmus,
     main="Elbow plot, ME setting 1",
     xlab = expression(delta),
     ylab = "No. nonzero coefficients",
     type="l",lty=1)
points(c(0.12),y=MonteCarlo_NCmean_gmus[c(7)],
       col="red",pch=20)
lines(delta_grid,MonteCarlo_NCmean_gmus0,
      lty=2)
points(c(0.12),y=MonteCarlo_NCmean_gmus0[c(7)],
       col="red",pch=20)

legend(0.3, 35, legend=c("CGMUS", "Equal ME CGMUS"),
       lty=c(1,2), cex=0.7)



##############

delta_esti_gmus<- 0.12#(0.14,0.2)   
delta_esti_gmus0<- 0.12#(0.14,0.2)   
delta_esti_gmus00<- 0.06 #0.06 or 0.08 from paper
delta_esti_gmus00_delta008<- 0.08 #0.06 or 0.08 from paper
## simulation, 100 dataset (each dataset use same set of coveriate-specific ME, sdX ) 
 
gmus_beta<-c()
gmus0_beta<-c()
gds_beta<-c()
lasso_beta<-c()
gmus00_beta<-c()
gmus00_beta_delta008<-c()
B=100
    
for (bstar in 1:B){
    
  set.seed(20211215+10*bstar)
  ll=generate_logistic_data(n=n,p=p,s=s,beta=beta,sdU=sdU,sdX=sdX)
    
  ## gmus
  merror=ll$sigmaUU
  gmus_fit2 <- gmus(ll$W, ll$y,merror=merror,
                    delta=delta_esti_gmus, family = "binomial")
  gmus_beta<-rbind(gmus_beta,t(gmus_fit2$beta))
  
  ## logistic lasso, lambda.min
  logistic_lasso_fit<-glmnet::glmnet(x=ll$W, y=ll$y, family = "binomial",
                                     lambda=gmus_fit2$lambda,alpha=1)
  lasso_beta<-rbind(lasso_beta,t(logistic_lasso_fit$beta))
  
  #gmus0
  merror_equal=diag(rep(sqrt(mean(diag(ll$sigmaUU)^2)),times=length(diag(ll$sigmaUU))))
  gmus0_fit2 <- gmus(ll$W, ll$y,merror=merror_equal, lambda= gmus_fit2$lambda,
                     delta=delta_esti_gmus0, family = "binomial")
  gmus0_beta<-rbind(gmus0_beta,t(gmus0_fit2$beta))
  
  #gds 
  gds_fit<-hdme::gds(ll$W,ll$y,lambda= gmus_fit2$lambda, family = "binomial")
  gds_beta<-rbind(gds_beta,t(gds_fit$beta))

  #gmus00 delta=0.06
  gmus00_fit<-hdme::gmus(ll$W,ll$y,delta=delta_esti_gmus00,
                         lambda = gmus_fit2$lambda,family="binomial")
  gmus00_beta<-rbind(gmus00_beta,t(gmus00_fit$beta))
  
  #gmus00 delta=0.08
  gmus00_fit<-hdme::gmus(ll$W,ll$y,delta=delta_esti_gmus00_delta008,
                         lambda = gmus_fit2$lambda,family="binomial")
  gmus00_beta_delta008<-rbind(gmus00_beta_delta008,t(gmus00_fit$beta))

}
  

#for (bstar in 1:B){ #gmus00 delta=0.06
#  
#  set.seed(20211215+10*bstar)
#  ll=generate_logistic_data(n=n,p=p,s=s,beta=beta,sdU=sdU,sdX=sdX)
#  gmus00_fit<-hdme::gmus(ll$W,ll$y,delta=delta_esti_gmus00,family="binomial")
#  gmus00_beta<-rbind(gmus00_beta,t(gmus00_fit$beta))
#  
#}

#for (bstar in 1:B){ #gmus00 delta=0.08
#  
#  set.seed(20211215+10*bstar)
#  ll=generate_logistic_data(n=n,p=p,s=s,beta=beta,sdU=sdU,sdX=sdX)
#  gmus00_fit<-hdme::gmus(ll$W,ll$y,delta=delta_esti_gmus00_delta008,family="binomial")
#  gmus00_beta_delta008<-rbind(gmus00_beta_delta008,t(gmus00_fit$beta))
#  
#}

  
#save(gmus_beta,file="gmus_beta_me02_equalME_lambdamin")
#save(gmus0_beta,file="gmus0_beta_me02_equalME_lambdamin")
#save(gmus00_beta,file="gmus00_beta_me02_equalME_lambdamin")
#save(gmus00_beta_delta008,file="gmus00_beta_delta008_me02_equalME_lambdamin")
#save(gds_beta,file="gds_beta_me02_equalME_lambdamin")
#save(lasso_beta,file="lasso_beta_me02_equalME_lambdamin")

######## second run till here

load(file="gmus_beta_me02_equalME_lambdamin")
load(file="gmus0_beta_me02_equalME_lambdamin")
load(file="gmus00_beta_me02_equalME_lambdamin")
load(file="gds_beta_me02_equalME_lambdamin")
load(file="lasso_beta_me02_equalME_lambdamin")
load(file="gmus00_beta_delta008_me02_equalME_lambdamin")


#summary_estbeta<-function(coef_mat){
#  temp=apply(coef_mat, MARGIN = 1,
#             FUN = function(x){
#               temp<-rep(0,times=p)
#               temp[abs(x)>1e-20]<- 1
#               true_beta_ind<-c(rep(1,s),rep(0,p-s))
#               
#               TP<- sum(temp[1:s]==true_beta_ind[1:s]) #true positive, number of coef correctly predicted !=0
#               TN<- sum(temp[(s+1):p]==true_beta_ind[(s+1):p]) # true negative, number of coef correctky predicted =0
#               FP<- sum(temp[(s+1):p]!=true_beta_ind[(s+1):p]) #false positive, predicted as nonzero coef, but actually zero
#               FN<- sum(temp[1:s]!=true_beta_ind[1:s]) #false negative
#               
#               sensitivity<-TP/(TP+FN)
#               specificity<-TN/(TN+FP)
#               accuracy<-(TP+TN)/(TP+TN+FP+FN)
#               precision<-TP/(TP+FP)
#               return(c(TP,TN,FP,FN,sensitivity,specificity,accuracy,precision))
#             }
#  )
#  output=as.data.frame(t(temp))
#  colnames(output)<-c("TP","TN","FP","FN","Sensitivity","Specificity","Accuracy","Precision")
#  return(output)
#}

n=200
p=500
s=10
beta <- c(rep(1,times=s),rep(0,times=p-s))


summary_estbeta<-function(coef_mat){
  temp=apply(coef_mat, MARGIN = 1,
             FUN = function(x){
               temp<-rep(0,times=p)
               temp[abs(x)>1e-20]<- 1
               true_beta_ind<-c(rep(1,s),rep(0,p-s))
               
               TP<- sum(temp[1:s]==true_beta_ind[1:s]) #true positive, number of coef correctly predicted !=0
               TN<- sum(temp[(s+1):p]==true_beta_ind[(s+1):p]) # true negative, number of coef correctky predicted =0
               FP<- sum(temp[(s+1):p]!=true_beta_ind[(s+1):p]) #false positive, predicted as nonzero coef, but actually zero
               FN<- sum(temp[1:s]!=true_beta_ind[1:s]) #false negative
               
               l1_error<-sum(abs(x-beta))
               l2_error<-sqrt(sum((x-beta)^2))
               
               sensitivity<-TP/(TP+FN)
               specificity<-TN/(TN+FP)
               accuracy<-(TP+TN)/(TP+TN+FP+FN)
               precision<-TP/(TP+FP)
               return(c(l1_error,l2_error,TP,TN,FP,FN,sensitivity,specificity,accuracy,precision))
             }
  )
  output=as.data.frame(t(temp))
  colnames(output)<-c("l1_error","l2_error","TP","TN","FP","FN","Sensitivity","Specificity","Accuracy","Precision")
  return(output)
}



#TP/(TP+FP) precision
mean(summary_estbeta(gmus_beta)[,8],na.rm=T)
mean(summary_estbeta(gmus0_beta)[,8],na.rm=T)
mean(summary_estbeta(gmus00_beta)[,8],na.rm=T)
mean(summary_estbeta(gmus00_beta_delta008)[,8],na.rm=T)
mean(summary_estbeta(gds_beta)[,8],na.rm=T)
mean(summary_estbeta(lasso_beta)[,8],na.rm=T)


colMeans(summary_estbeta(gmus_beta)[summary_estbeta(gmus_beta)[,3]!=0,])
colMeans(summary_estbeta(gmus0_beta)[summary_estbeta(gmus0_beta)[,3]!=0,])
colMeans(summary_estbeta(gmus00_beta)[summary_estbeta(gmus00_beta)[,3]!=0,])
colMeans(summary_estbeta(gmus00_beta_delta008)[summary_estbeta(gmus00_beta_delta008)[,3]!=0,])
colMeans(summary_estbeta(gds_beta)[summary_estbeta(gds_beta)[,3]!=0,])
colMeans(summary_estbeta(lasso_beta)[summary_estbeta(lasso_beta)[,3]!=0,])































#mean TP
mean(summary_estbeta(gmus_beta)[,1],na.rm=T)
mean(summary_estbeta(gmus0_beta)[,1],na.rm=T)
mean(summary_estbeta(gmus00_beta)[,1],na.rm=T)
mean(summary_estbeta(gmus00_beta_delta008)[,1],na.rm=T)
mean(summary_estbeta(gds_beta)[,1],na.rm=T)
mean(summary_estbeta(lasso_beta)[,1],na.rm=T)

#mean TN
mean(summary_estbeta(gmus_beta)[,2],na.rm=T)
mean(summary_estbeta(gmus0_beta)[,2],na.rm=T)
mean(summary_estbeta(gmus00_beta)[,2],na.rm=T)
mean(summary_estbeta(gmus00_beta_delta008)[,2],na.rm=T)
mean(summary_estbeta(gds_beta)[,2],na.rm=T)
mean(summary_estbeta(lasso_beta)[,2],na.rm=T)

#mean Fp
mean(summary_estbeta(gmus_beta)[,3],na.rm=T)
mean(summary_estbeta(gmus0_beta)[,3],na.rm=T)
mean(summary_estbeta(gmus00_beta)[,3],na.rm=T)
mean(summary_estbeta(gmus00_beta_delta008)[,3],na.rm=T)
mean(summary_estbeta(gds_beta)[,3],na.rm=T)
mean(summary_estbeta(lasso_beta)[,3],na.rm=T)


#mean sensitivity
mean(summary_estbeta(gmus_beta)[,5],na.rm=T)
mean(summary_estbeta(gmus0_beta)[,5],na.rm=T)
mean(summary_estbeta(gmus00_beta)[,5],na.rm=T)
mean(summary_estbeta(gmus00_beta_delta008)[,5],na.rm=T)
mean(summary_estbeta(gds_beta)[,5],na.rm=T)
mean(summary_estbeta(lasso_beta)[,5],na.rm=T)
