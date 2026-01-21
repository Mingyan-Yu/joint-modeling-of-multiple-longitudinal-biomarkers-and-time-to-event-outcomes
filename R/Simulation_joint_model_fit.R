###### Joint modeling of longitudinal biomarkers and time-to-event outcomes ######

## load packages
library(nimble)
library(data.table)
library(coda)
library(parallel)
library(nlme)
library(GLMMadaptive)
library(dplyr)

k <- Sys.getenv("SLURM_ARRAY_TASK_ID")
print(k)

set.seed(20251210+as.numeric(k)*100)

## create 3 parallel clusters
cl <- makeCluster(3)

## load simulated data
load(paste0(k, ".result/LOD.Rdata"))
load(paste0(k, ".result/surv_biomarker_data.Rdata"))
load(paste0(k, ".result/t.Rdata")
load(paste0(k, ".result/n.Rdata")

## a user-defined function to combine the Cholesky of the correlation matrix
## with a vector of standard deviations to produce the Cholesky of the covariance matrix
uppertri_mult_diag <- nimbleFunction(
  run = function(mat = double(2), vec = double(1)) {
    returnType(double(2))
    p <- length(vec)
    out <- matrix(nrow = p, ncol = p, init = FALSE)
    for(i in 1:p)
      out[ , i] <- mat[ , i] * vec[i]
    return(out)
  })

## create user-defined Inverse-Gaussian distribution functions
## density function
dinversegauss <- nimbleFunction(
  run = function(x = double(0), mean_mu = double(0, default = 1),
                 lambda = double(0, default = 1),
                 log = integer(0, default = 0)){
    returnType(double(0))
    logProb <- 0.5*(log(lambda) - log(2*pi*x^3)) - lambda*(x-mean_mu)^2/(2*mean_mu^2*x)
    if(log) return(logProb)
    else return(exp(logProb))
  }
)

## random number generation function
rinversegauss <- nimbleFunction(
  run = function(n = integer(0), mean_mu = double(0, default = 1),
                 lambda = double(0, default = 1)){
    returnType(double(0))
    y <- rchisq(1,1)
    x <- mean_mu + (mean_mu^2*y)/(2*lambda) - mean_mu/(2*lambda)*sqrt(4*mean_mu*lambda*y+(mean_mu*y)^2)
    if(runif(1,0,1) <= mean_mu/(mean_mu+x)){
      return(x)
    } else {
      return(mean_mu^2/x)
    }
  }
)

## export the user-defined functions to the clusters
clusterExport(cl, "uppertri_mult_diag")
clusterExport(cl, c("dinversegauss", "rinversegauss"))
clusterExport(cl, "registerDistributions")
clusterEvalQ(cl, {
  library(nimble)
})

clusterEvalQ(cl, {
  registerDistributions(list(
    dinversegauss = list(
      BUGSdist = "dinversegauss(mean_mu, lambda)",
      types = c("value=double(0)", "mean_mu=double(0)", "lambda=double(0)")
    )
  ))
})

## construct the NIMBLE joint model for two biomarkers with one subject to detection limits
## and time-to-event outcomes
## with time on the original scale
joint_model <- nimbleCode({
  for(i in 1:N){
    for(j in 1:n[i]){
      ## specify mean functions as a combination of fixed effects and random effects
      mean_marker[i,j,1] <- beta1[1]+beta1[2]*time[i,j]+
        beta1[3]*time[i,j]^2+beta1[4]*time[i,j]^3+
        b1[i,1]+b1[i,2]*time[i,j]
      
      mean_marker[i,j,2] <- beta2[1]+beta2[2]*time[i,j]+
        beta2[3]*time[i,j]^2+b2[i,1]+b2[i,2]*time[i,j]
    }
  }
  
  ## model biomarker values that were either observed or censored
  for(i in 1:N){
    for(j in 1:n[i]){
      ## model biomarker1 (fully observed) by its marginal normal distribution
      marker1[i,j] ~ dnorm(mean_marker[i,j,1], sd = sqrt(var_sq1[i]))

      ## model biomarker2 (with left-censoring) by a conditional normal distribution
      ## conditioned on the value of biomarker1
      ## with imputation incorporated for values below LOD
      censored_marker[i,j] ~ dinterval(marker2[i,j], LOD)
      marker2[i,j] ~ dnorm(mean_marker[i,j,2]+r[i]*sqrt(var_sq2[i]/var_sq1[i])*
                             (marker1[i,j]-mean_marker[i,j,1]),
                           sd = sqrt(1-r[i]^2)*sqrt(var_sq2[i]))
    }
  }
  
  ## priors for fixed effects
  for(i in 1:4){
    beta1[i] ~ dnorm(0, sd = 10)
  }
  for(i in 1:3){
    beta2[i] ~ dnorm(0, sd = 10)
  }
  
  ## priors for random effects
  ### lkj priors for the cholesky decomposition of the correlation matrix between random effects
  Ucorr_b1[1:2, 1:2] ~ dlkj_corr_cholesky(1.0, 2)  ## biomarker1
  Ucorr_b2[1:2, 1:2] ~ dlkj_corr_cholesky(1.0, 2)  ## biomarker2
  for(i in 1:2){
    ### Half-Cauchy prior for the standard deviations of random effects
    sigma_b1[i] ~ T(dt(0, pow(2.5,-2), 1),0,)  ## biomarker1
    sigma_b2[i] ~ T(dt(0, pow(2.5,-2), 1),0,)  ## biomarker2
  }
  ### obtain the cholesky decomposition of the covariance matrix between random effects
  Ucov_b1[1:2, 1:2] <- uppertri_mult_diag(Ucorr_b1[1:2,1:2], sigma_b1[1:2])  ## biomarker1
  Ucov_b2[1:2, 1:2] <- uppertri_mult_diag(Ucorr_b2[1:2,1:2], sigma_b2[1:2])  ## biomarker2
  for(i in 1:2){
    mean_b_prior[i] <- 0
  }
  for(i in 1:N){
    ## multivariate-Gaussian priors for the random effects
    b1[i,1:2] ~ dmnorm(mean_b_prior[1:2], cholesky = Ucov_b1[1:2,1:2], prec_param = 0)
    
    b2[i,1:2] ~ dmnorm(mean_b_prior[1:2], cholesky = Ucov_b2[1:2,1:2], prec_param = 0)
  }
  
  ## priors for residual variance
  gamma1 ~ dnorm(0, sd = 10)
  Psi1 ~ T(dt(0, pow(2.5,-2), 1),0,)
  gamma2 ~ dnorm(0, sd = 10)
  Psi2 ~ T(dt(0, pow(2.5,-2), 1),0,)
  for(i in 1:N){
    var_sq1[i] ~ dlnorm(gamma1, sd = Psi1)
    var_sq2[i] ~ dlnorm(gamma2, sd = Psi2)
  }
  
  ## priors for covariance matrix
  for(i in 1:N){
    t_r[i] ~ dbeta(a, b)
    r[i] <- t_r[i]*2-1
  }
  a ~ dexp(0.1)
  b ~ dexp(0.1)
  
  ## create the matrix of biomarker predictors
  for(i in 1:N){
    biomarker_predictor[i,1] <- 1
    biomarker_predictor[i,2] <- (b1[i,1]-0)/sigma_b1[1]
    biomarker_predictor[i,3] <- (b1[i,2]-0)/sigma_b1[2]
    biomarker_predictor[i,4] <- (log(var_sq1[i])-gamma1)/Psi1
    biomarker_predictor[i,5] <- (b2[i,1]-0)/sigma_b2[1]
    biomarker_predictor[i,6] <- (b2[i,2]-0)/sigma_b2[2]
    biomarker_predictor[i,7] <- (log(var_sq2[i])-gamma2)/Psi2
  }
  
  ## start to model time-to-event outcomes
  for(i in 1:N){
    mean_mu[i] <- -exp(inprod(biomarker_predictor[i,1:7],alpha[1:7]))/inprod(biomarker_predictor[i,1:7],eta[1:7])
    lambda[i] <- exp(inprod(biomarker_predictor[i,1:7],alpha[1:7]))^2
    
    censored_event[i] ~ dinterval(t[i], c[i])
    t[i] ~ dinversegauss(mean_mu[i], lambda[i])
  }
  
  for(j in 1:7){
    alpha[j] ~ dnorm(0, sd = 10)
    eta[j] ~ dnorm(0, sd = 10)
  }
})


















     

