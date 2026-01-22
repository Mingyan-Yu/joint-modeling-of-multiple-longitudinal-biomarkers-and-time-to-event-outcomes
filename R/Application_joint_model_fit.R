###### Apply the joint model to the application dataset ######

## load packages
library(nimble)
library(coda)
library(parallel)
library(nlme)
library(GLMMadaptive)
library(dplyr)

k <- Sys.getenv("SLURM_ARRAY_TASK_ID")

set.seed(20260115+as.numeric(k)*2026)

## create 3 parallel clusters
cl <- makeCluster(3)

## import the application dataset
load("swan_restrict4_both.Rdata")

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

## export the user-defined function to the clusters
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
## with measurement time standardized
joint_nimble <- nimbleCode({
  for(i in 1:N){
    for(j in 1:n[i]){
      ## specify mean functions as a combination of fixed effects and random effects
      ## with standardized measurement times
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
  
  ## transform beta1 back to correspond to the original time scale
  beta1_ori[1] <- beta1[1]-beta1[2]/sd_t*mean_t+beta1[3]/sd_t^2*mean_t^2-beta1[4]/sd_t^3*mean_t^3
  beta1_ori[2] <- beta1[2]/sd_t-2*beta1[3]/sd_t^2*mean_t+3*beta1[4]/sd_t^3*mean_t^2
  beta1_ori[3] <- beta1[3]/sd_t^2-3*beta1[4]/sd_t^3*mean_t
  beta1_ori[4] <- beta1[4]/sd_t^3
  
  ## transform beta2 back to correspond to the original time scale
  beta2_ori[1] <- beta2[1]-beta2[2]/sd_t*mean_t+beta2[3]/sd_t^2*mean_t^2
  beta2_ori[2] <- beta2[2]/sd_t-2*beta2[3]/sd_t^2*mean_t
  beta2_ori[3] <- beta2[3]/sd_t^2
  
  ## priors for random effects
  ### lkj priors for the cholesky decomposition of the correlation matrix between random effects
  Ucorr_b1[1:2, 1:2] ~ dlkj_corr_cholesky(1.0, 2)  ## biomarker1
  Ucorr_b2[1:2, 1:2] ~ dlkj_corr_cholesky(1.0, 2)  ## biomarker2
  ### calculate the correlation between random effects
  Omega_b1 <- inprod(Ucorr_b1[1:2, 1], Ucorr_b1[1:2, 2])  ## biomarker1
  Omega_b2 <- inprod(Ucorr_b2[1:2, 1], Ucorr_b2[1:2, 2])  ## biomarker2
  for(i in 1:2){
    ### Half-Cauchy prior for the standard deviations of random effects
    sigma_b1[i] ~ T(dt(0, pow(2.5,-2), 1),0,)
    sigma_b2[i] ~ T(dt(0, pow(2.5,-2), 1),0,)
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
    
    ## transform b1 and b2 back to correspond to the original time scale
    b1_ori[i,1] <- b1[i,1] - b1[i,2]*mean_t/sd_t
    b1_ori[i,2] <- b1[i,2]/sd_t
    b2_ori[i,1] <- b2[i,1] - b2[i,2]*mean_t/sd_t
    b2_ori[i,2] <- b2[i,2]/sd_t
  }
  
  ## transform sigma_b1 and sigma_b2 back to correspond to the original time scale
  sigma_b1_ori[1] <- sqrt(sigma_b1[1]^2-2*mean_t/sd_t*sigma_b1[1]*sigma_b1[2]*Omega_b1+
                            mean_t^2/sd_t^2*sigma_b1[2]^2)
  sigma_b1_ori[2] <- sigma_b1[2]/sd_t
  sigma_b2_ori[1] <- sqrt(sigma_b2[1]^2-2*mean_t/sd_t*sigma_b2[1]*sigma_b2[2]*Omega_b2+
                            mean_t^2/sd_t^2*sigma_b2[2]^2)
  sigma_b2_ori[2] <- sigma_b2[2]/sd_t
  
  ## priors for residual variance
  gamma1 ~ dnorm(0, sd = 10)
  Psi1 ~ T(dt(0, pow(2.5,-2), 1),0,)
  gamma2 ~ dnorm(0, sd = 10)
  Psi2 ~ T(dt(0, pow(2.5,-2), 1),0,)
  for(i in 1:N){
    var_sq1[i] ~ dlnorm(gamma1, sd = Psi1)
    var_sq2[i] ~ dlnorm(gamma2, sd = Psi2)
  }
  
  ## priors for residual correlation
  for(i in 1:N){
    t_r[i] ~ dbeta(a, b)
    r[i] <- t_r[i]*2-1
  }
  a ~ dexp(0.1)
  b ~ dexp(0.1)
  
  ## create the matrix of biomarker predictors with standardization
  for(i in 1:N){
    biomarker_predictor[i,1] <- 1
    biomarker_predictor[i,2] <- (b1_ori[i,1]-0)/sigma_b1_ori[1]
    biomarker_predictor[i,3] <- (b1_ori[i,2]-0)/sigma_b1_ori[2]
    biomarker_predictor[i,4] <- (log(var_sq1[i])-gamma1)/Psi1
    biomarker_predictor[i,5] <- (b2_ori[i,1]-0)/sigma_b2_ori[1]
    biomarker_predictor[i,6] <- (b2_ori[i,2]-0)/sigma_b2_ori[2]
    biomarker_predictor[i,7] <- (log(var_sq2[i])-gamma2)/Psi2
    biomarker_predictor[i,8] <- BMI[i]
  }
  
  ## start to model time-to-event outcomes
  for(i in 1:N){
    mean_mu[i] <- -exp(inprod(biomarker_predictor[i,1:8],alpha[1:8]))/inprod(biomarker_predictor[i,1:8],eta[1:8])
    lambda[i] <- exp(inprod(biomarker_predictor[i,1:8],alpha[1:8]))^2
    
    censored_event[i] ~ dinterval(t[i], c[i])
    t[i] ~ dinversegauss(mean_mu[i], lambda[i])
  }
  
  for(j in 1:8){
    alpha[j] ~ dnorm(0, sd = 10)
    eta[j] ~ dnorm(0, sd = 10)
  }
})

LOD <- log(1.85)
## construct data for fitting the joint model
## turn left-censored AMH values (log1.85 or log1.45) into NA
marker2 <- ifelse(swanDat_restrict4$marker2 == log(1.85) |
                    swanDat_restrict4$marker2 == log(1.45), NA, 
                  swanDat_restrict4$marker2)
## construct a censoring indicator for AMH
censored_marker <- marker2
censored_marker[] <- as.numeric(!is.na(marker2))

joint_data <- list(marker1 = swanDat_restrict4$marker1,
                   marker2 = marker2,
                   censored_marker = censored_marker,
                   censored_event = as.numeric(swanDat_restrict4$swan_restrict4_unique$status==0),
                   t = ifelse(swanDat_restrict4$swan_restrict4_unique$status==1,
                              swanDat_restrict4$swan_restrict4_unique$event_time-42, NA)) ## subtract event times by 42 to align with biomarker measurement times

## standardize measurement times
mean_t = mean(swanDat_restrict4$t[swanDat_restrict4$t!=0])
sd_t = sd(swanDat_restrict4$t[swanDat_restrict4$t!=0])
t_std = (swanDat_restrict4$t-mean_t)/sd_t

## construct constants for fitting the joint model
joint_constants <- list(N = swanDat_restrict4$N, n = swanDat_restrict4$n, 
                        time = t_std, LOD = LOD, mean_t = mean_t, sd_t = sd_t,
                        c = ifelse(swanDat_restrict4$swan_restrict4_unique$status==1, Inf,
                                   swanDat_restrict4$swan_restrict4_unique$event_time-42),
                        BMI = swanDat_restrict4$BMI)

## specify initial values for some parameters
joint_inits <- list()

## for beta1, sigma_b1, Omega_b1, gamma1 and Psi1: 
## fit a linear mixed effects model on FSH
FSH_data <- data.frame()
for(i in 1:swanDat_restrict4$N){
  FSH_i <- data.frame("ID" = rep(i, swanDat_restrict4$n[i]),
                      "Time" = t_std[i,1:swanDat_restrict4$n[i]],
                      "FSH_log" = swanDat_restrict4$marker1[i,1:swanDat_restrict4$n[i]])
  FSH_data <- rbind(FSH_data, FSH_i)
}

FSH_fit <- lme(fixed = FSH_log ~ Time + I(Time^2) + I(Time^3),
               random = ~ 1 + Time | ID,
               data = FSH_data)

joint_inits$beta1 <- fixef(FSH_fit) ## from fixed effects
joint_inits$sigma_b1 <- apply(ranef(FSH_fit), MARGIN = 2, FUN = sd) ## from sd of random effects
joint_inits$Ucorr_b1 <- rlkj_corr_cholesky(n = 1, eta = 1, p = 2) ## randomly generate from a lkj distribution
joint_inits$Ucov_b1 <- uppertri_mult_diag(joint_inits$Ucorr_b1[1:2,1:2], 
                                          joint_inits$sigma_b1[1:2])
FSH_data$residual <- FSH_data$FSH_log - FSH_fit$fitted[,2]  ## obtain individual-specific residuals
## calculate individual-specific residual variance
var_FSH <- FSH_data %>%
  group_by(ID) %>%
  summarise(Si1 = var(residual))
joint_inits$gamma1 <- mean(log(var_FSH$Si1)) ## from mean of log-transformed residual variance
joint_inits$Psi1 <- sd(log(var_FSH$Si1)) #from sd of log-transformed residual variance

joint_inits$var_sq1 <- rlnorm(joint_constants$N, meanlog = joint_inits$gamma1,
                                      sdlog = joint_inits$Psi1)

## for beta2, sigma_b2, Omega_b2, gamma2 and Psi2:
## fit a tobit mixed effects model on AMH
AMH_data <- data.frame()
for(i in 1:swanDat_restrict4$N){
  AMH_i <- data.frame("ID" = rep(i, swanDat_restrict4$n[i]),
                      "Time" = t_std[i,1:swanDat_restrict4$n[i]],
                      "AMH_log" = swanDat_restrict4$marker2[i,1:swanDat_restrict4$n[i]])
  AMH_data <- rbind(AMH_data, AMH_i)
}
AMH_data$AMH_log <- ifelse(AMH_data$AMH_log==log(1.45), LOD, AMH_data$AMH_log)

AMH_data$censor <- ifelse(AMH_data$AMH_log == LOD, 1, 0)

AMH_fit <- mixed_model(fixed = cbind(AMH_log, censor) ~ Time + I(Time^2),
                       random = ~ 1 + Time | ID,
                       data = AMH_data,
                       family = censored.normal())

joint_inits$beta2 <- fixef(AMH_fit) ## from fixed effects
joint_inits$sigma_b2 <- apply(ranef(AMH_fit), MARGIN = 2, FUN = sd) ## from sd of random effects
joint_inits$Ucorr_b2 <- rlkj_corr_cholesky(n = 1, eta = 1, p = 2) ## randomly generate from a lkj distribution
joint_inits$Ucov_b2 <- uppertri_mult_diag(joint_inits$Ucorr_b2[1:2,1:2], 
                                          joint_inits$sigma_b2[1:2])
## obtain individual-specific residuals
AMH_data$residual <- AMH_data$AMH_log - GLMMadaptive::fitted(AMH_fit, type = "subject_specific")
## calculate individual-specific residual variance
var_AMH <- AMH_data %>%
  group_by(ID) %>%
  summarise(Si2 = var(residual))
joint_inits$gamma2 <- mean(log(var_AMH$Si2)) ## from mean of log-transformed residual variance
joint_inits$Psi2 <- sd(log(var_AMH$Si2)) #from sd of log-transformed residual variance

joint_inits$var_sq2 <- rlnorm(joint_constants$N, meanlog = joint_inits$gamma2,
                              sdlog = joint_inits$Psi2)

## for a and b: calculate the correlation between FSH residuals and AMH residuals
marker_cor_data <- merge(FSH_data, AMH_data,
                         by.x = c("ID", "Time"),
                         by.y = c("ID", "Time"),
                         all.y = TRUE)
marker_cor <- marker_cor_data %>%
  group_by(ID) %>%
  summarise(marker_cor = cor(residual.x, residual.y))
marker_cor_transform <- (marker_cor$marker_cor+1)/2

joint_inits$a <- mean(marker_cor_transform, na.rm = TRUE)^2*
  (1-mean(marker_cor_transform, na.rm = TRUE))/var(marker_cor_transform, na.rm = TRUE) - 
  mean(marker_cor_transform, na.rm = TRUE)
joint_inits$b <- joint_inits$a*(1/mean(marker_cor_transform, na.rm = TRUE)-1)

print(joint_inits$a)
print(joint_inits$b)

joint_inits$t_r <- rbeta(joint_constants$N, joint_inits$a, joint_inits$b)
joint_inits$r <- joint_inits$t_r*2-1

## specify the rest of the initial values
## for AMH below LOD, specify the initial values for the imputation as log(1.45)
joint_inits$marker2 <- ifelse(is.na(marker2), log(1.45), NA)

joint_inits$alpha <- rnorm(n = 8, mean = 0, sd = 0.5)
joint_inits$eta <- rnorm(n = 8, mean = 0, sd = 0.5)
## for right-censored age at FMP, specify the initial values for imputation as censoring time + 0.3
joint_inits$t <- rep(NA, joint_constants$N)
joint_inits$t[joint_data$censored_event==1] <- joint_constants$c[joint_data$censored_event==1] + 0.3

## start to build and fit the NIMBLE model in parallel
run_MCMC_parallel <- function(seed, data, constants, code, inits){
  library(nimble)
  
  ## build the nimble model
  joint_model <- nimbleModel(code = code, name = "joint_model",
                             data = data, constants = constants,
                             inits = inits)
  
  ## configure the model
  joint_model.config <- configureMCMC(joint_model, print = TRUE)
  joint_model.config$addMonitors("b1_ori", "b2_ori", "Omega_b1", "Omega_b2",
                                 "beta1_ori", "beta2_ori",
                                 "sigma_b1_ori", "sigma_b2_ori")
  
  ## change samplers for some parameters
  ### replace samplers for beta1 and beta2 with automated factor slice samplers
  joint_model.config$removeSamplers("beta1[]", "beta2[]")
  joint_model.config$addSampler(target = c("beta1[]", "beta2[]"),
                                type = "AF_slice")
  
  ## replace samplers for alpha and eta with automated factor slice sampler
  joint_model.config$removeSamplers("alpha[]", "eta[]")
  joint_model.config$addSampler(target = c("alpha[]", "eta[]"),
                                type = "AF_slice")
  
  ### replace samplers for some other parameters with slice samplers
  joint_model.config$removeSamplers(c("sigma_b1[1]", "sigma_b2[1]",
                                      "sigma_b1[2]", "sigma_b2[2]",
                                      "Psi1", "Psi2", "a", "b"))
  joint_model.config$addSampler(target = "sigma_b1[]", type = "AF_slice")
  joint_model.config$addSampler(target = "sigma_b2[]", type = "AF_slice")
  joint_model.config$addSampler(target = "Psi1",
                                type = "slice")
  joint_model.config$addSampler(target = "Psi2",
                                type = "slice")
  joint_model.config$addSampler(target = "a",
                                type = "slice")
  joint_model.config$addSampler(target = "b",
                                type = "slice")
  
  ## build an MCMC object for this model
  joint_model.MCMC <- buildMCMC(joint_model.config)
  
  ## compile the model and the MCMC object
  Com_joint.Model <- compileNimble(joint_model)
  Com_joint.MCMC <- compileNimble(joint_model.MCMC, project = joint_model)
  
  Sys.time()
  ## run the model
  run_joint_samples <- runMCMC(Com_joint.MCMC, niter = 240000,
                               nburnin = 120000, progressBar = TRUE,
                               summary = TRUE, setSeed = seed)
  Sys.time()
  
  return(run_joint_samples)
}

chain_output <- parLapply(cl = cl, X = 20260130+as.numeric(k)+(1:3)*100,
                          fun = run_MCMC_parallel,
                          data = joint_data,
                          code = joint_nimble,
                          constants = joint_constants,
                          inits = joint_inits)

save(chain_output, file = "chain_output.Rdata")

col <- colnames(chain_output[[1]]$samples)
n.sims <- nrow(chain_output[[1]]$samples)

## post-sampling processing
## calculate correlation matrix for the random effects
Omega_b1_sample <- list()
for(i in 1:3){
  cov_b1_sample <- chain_output[[i]]$samples[,"Omega_b1"] * 
    chain_output[[i]]$samples[,"sigma_b1[1]"] * chain_output[[i]]$samples[,"sigma_b1[2]"]
  cov_b1_ori_sample <- cov_b1_sample/sd_t - mean_t/sd_t^2*chain_output[[i]]$samples[,"sigma_b1[2]"]^2
  Omega_b1_sample[[i]] <- cov_b1_ori_sample/chain_output[[i]]$samples[,"sigma_b1_ori[1]"]/
    chain_output[[i]]$samples[,"sigma_b1_ori[2]"]
}

Omega_b2_sample <- list()
for(i in 1:3){
  cov_b2_sample <- chain_output[[i]]$samples[,"Omega_b2"] * 
    chain_output[[i]]$samples[,"sigma_b2[1]"] * chain_output[[i]]$samples[,"sigma_b2[2]"]
  cov_b2_ori_sample <- cov_b2_sample/sd_t - mean_t/sd_t^2*chain_output[[i]]$samples[,"sigma_b2[2]"]^2
  Omega_b2_sample[[i]] <- cov_b2_ori_sample/chain_output[[i]]$samples[,"sigma_b2_ori[1]"]/
    chain_output[[i]]$samples[,"sigma_b2_ori[2]"]
}

## make a traceplot of the posterior samples
pdf(paste0(k, ".result/traceplot.pdf"))
coda::traceplot(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^beta1_ori",col)]),
                          mcmc(chain_output[[2]]$samples[,grep("^beta1_ori",col)]),
                          mcmc(chain_output[[3]]$samples[,grep("^beta1_ori",col)])))
coda::traceplot(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^Psi1",col)]),
                          mcmc(chain_output[[2]]$samples[,grep("^Psi1",col)]),
                          mcmc(chain_output[[3]]$samples[,grep("^Psi1",col)])))
coda::traceplot(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^gamma1",col)]),
                          mcmc(chain_output[[2]]$samples[,grep("^gamma1",col)]),
                          mcmc(chain_output[[3]]$samples[,grep("^gamma1",col)])))
coda::traceplot(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^sigma_b1_ori",col)]),
                          mcmc(chain_output[[2]]$samples[,grep("^sigma_b1_ori",col)]),
                          mcmc(chain_output[[3]]$samples[,grep("^sigma_b1_ori",col)])))
coda::traceplot(mcmc.list(mcmc(Omega_b1_sample[[1]]),
                          mcmc(Omega_b1_sample[[2]]),
                          mcmc(Omega_b1_sample[[3]])))
coda::traceplot(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^beta2_ori",col)]),
                          mcmc(chain_output[[2]]$samples[,grep("^beta2_ori",col)]),
                          mcmc(chain_output[[3]]$samples[,grep("^beta2_ori",col)])))
coda::traceplot(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^Psi2",col)]),
                          mcmc(chain_output[[2]]$samples[,grep("^Psi2",col)]),
                          mcmc(chain_output[[3]]$samples[,grep("^Psi2",col)])))
coda::traceplot(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^gamma2",col)]),
                          mcmc(chain_output[[2]]$samples[,grep("^gamma2",col)]),
                          mcmc(chain_output[[3]]$samples[,grep("^gamma2",col)])))
coda::traceplot(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^sigma_b2_ori",col)]),
                          mcmc(chain_output[[2]]$samples[,grep("^sigma_b2_ori",col)]),
                          mcmc(chain_output[[3]]$samples[,grep("^sigma_b2_ori",col)])))
coda::traceplot(mcmc.list(mcmc(Omega_b2_sample[[1]]),
                          mcmc(Omega_b2_sample[[2]]),
                          mcmc(Omega_b2_sample[[3]])))
coda::traceplot(mcmc.list(mcmc(chain_output[[1]]$samples[,"a"]),
                          mcmc(chain_output[[2]]$samples[,"a"]),
                          mcmc(chain_output[[3]]$samples[,"a"])))
coda::traceplot(mcmc.list(mcmc(chain_output[[1]]$samples[,"b"]),
                          mcmc(chain_output[[2]]$samples[,"b"]),
                          mcmc(chain_output[[3]]$samples[,"b"])))

coda::traceplot(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^alpha",col)]),
                          mcmc(chain_output[[2]]$samples[,grep("^alpha",col)]),
                          mcmc(chain_output[[3]]$samples[,grep("^alpha",col)])))
coda::traceplot(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^eta",col)]),
                          mcmc(chain_output[[2]]$samples[,grep("^eta",col)]),
                          mcmc(chain_output[[3]]$samples[,grep("^eta",col)])))

dev.off()

## calculate posterior mean and 95% CIs
beta1_est <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^beta1_ori",col)]),
                               mcmc(chain_output[[2]]$samples[,grep("^beta1_ori",col)]),
                               mcmc(chain_output[[3]]$samples[,grep("^beta1_ori",col)])))$statistics[,1]
Psi1_est <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^Psi1",col)]),
                              mcmc(chain_output[[2]]$samples[,grep("^Psi1",col)]),
                              mcmc(chain_output[[3]]$samples[,grep("^Psi1",col)])))$statistics[1]
gamma1_est <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^gamma1",col)]),
                                mcmc(chain_output[[2]]$samples[,grep("^gamma1",col)]),
                                mcmc(chain_output[[3]]$samples[,grep("^gamma1",col)])))$statistics[1]
sigma_b1_est <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^sigma_b1_ori",col)]),
                                  mcmc(chain_output[[2]]$samples[,grep("^sigma_b1_ori",col)]),
                                  mcmc(chain_output[[3]]$samples[,grep("^sigma_b1_ori",col)])))$statistics[,1]
Omega_b1_est <- summary(mcmc.list(mcmc(Omega_b1_sample[[1]]),
                                  mcmc(Omega_b1_sample[[2]]),
                                  mcmc(Omega_b1_sample[[3]])))$statistics[1]
beta2_est <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^beta2_ori",col)]),
                               mcmc(chain_output[[2]]$samples[,grep("^beta2_ori",col)]),
                               mcmc(chain_output[[3]]$samples[,grep("^beta2_ori",col)])))$statistics[,1]
Psi2_est <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^Psi2",col)]),
                              mcmc(chain_output[[2]]$samples[,grep("^Psi2",col)]),
                              mcmc(chain_output[[3]]$samples[,grep("^Psi2",col)])))$statistics[1]
gamma2_est <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^gamma2",col)]),
                                mcmc(chain_output[[2]]$samples[,grep("^gamma2",col)]),
                                mcmc(chain_output[[3]]$samples[,grep("^gamma2",col)])))$statistics[1]
sigma_b2_est <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^sigma_b2_ori",col)]),
                                  mcmc(chain_output[[2]]$samples[,grep("^sigma_b2_ori",col)]),
                                  mcmc(chain_output[[3]]$samples[,grep("^sigma_b2_ori",col)])))$statistics[,1]
Omega_b2_est <- summary(mcmc.list(mcmc(Omega_b2_sample[[1]]),
                                  mcmc(Omega_b2_sample[[2]]),
                                  mcmc(Omega_b2_sample[[3]])))$statistics[1]
a_est <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,"a"]),
                           mcmc(chain_output[[2]]$samples[,"a"]),
                           mcmc(chain_output[[3]]$samples[,"a"])))$statistics[1]
b_est <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,"b"]),
                           mcmc(chain_output[[2]]$samples[,"b"]),
                           mcmc(chain_output[[3]]$samples[,"b"])))$statistics[1]
alpha_est <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^alpha",col)]),
                               mcmc(chain_output[[2]]$samples[,grep("^alpha",col)]),
                               mcmc(chain_output[[3]]$samples[,grep("^alpha",col)])))$statistics[,1]
eta_est <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^eta",col)]),
                             mcmc(chain_output[[2]]$samples[,grep("^eta",col)]),
                             mcmc(chain_output[[3]]$samples[,grep("^eta",col)])))$statistics[,1]


beta1_CI <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^beta1_ori",col)]),
                              mcmc(chain_output[[2]]$samples[,grep("^beta1_ori",col)]),
                              mcmc(chain_output[[3]]$samples[,grep("^beta1_ori",col)])))$quantiles[,c(1,5)]
Psi1_CI <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^Psi1",col)]),
                             mcmc(chain_output[[2]]$samples[,grep("^Psi1",col)]),
                             mcmc(chain_output[[3]]$samples[,grep("^Psi1",col)])))$quantiles[c(1,5)]
gamma1_CI <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^gamma1",col)]),
                               mcmc(chain_output[[2]]$samples[,grep("^gamma1",col)]),
                               mcmc(chain_output[[3]]$samples[,grep("^gamma1",col)])))$quantiles[c(1,5)]
sigma_b1_CI <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^sigma_b1_ori",col)]),
                                 mcmc(chain_output[[2]]$samples[,grep("^sigma_b1_ori",col)]),
                                 mcmc(chain_output[[3]]$samples[,grep("^sigma_b1_ori",col)])))$quantiles[,c(1,5)]
Omega_b1_CI <- summary(mcmc.list(mcmc(Omega_b1_sample[[1]]),
                                 mcmc(Omega_b1_sample[[2]]),
                                 mcmc(Omega_b1_sample[[3]])))$quantiles[c(1,5)]
beta2_CI <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^beta2_ori",col)]),
                              mcmc(chain_output[[2]]$samples[,grep("^beta2_ori",col)]),
                              mcmc(chain_output[[3]]$samples[,grep("^beta2_ori",col)])))$quantiles[,c(1,5)]
Psi2_CI <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^Psi2",col)]),
                             mcmc(chain_output[[2]]$samples[,grep("^Psi2",col)]),
                             mcmc(chain_output[[3]]$samples[,grep("^Psi2",col)])))$quantiles[c(1,5)]
gamma2_CI <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^gamma2",col)]),
                               mcmc(chain_output[[2]]$samples[,grep("^gamma2",col)]),
                               mcmc(chain_output[[3]]$samples[,grep("^gamma2",col)])))$quantiles[c(1,5)]
sigma_b2_CI <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^sigma_b2_ori",col)]),
                                 mcmc(chain_output[[2]]$samples[,grep("^sigma_b2_ori",col)]),
                                 mcmc(chain_output[[3]]$samples[,grep("^sigma_b2_ori",col)])))$quantiles[,c(1,5)]
Omega_b2_CI <- summary(mcmc.list(mcmc(Omega_b2_sample[[1]]),
                                 mcmc(Omega_b2_sample[[2]]),
                                 mcmc(Omega_b2_sample[[3]])))$quantiles[c(1,5)]
a_CI <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,"a"]),
                          mcmc(chain_output[[2]]$samples[,"a"]),
                          mcmc(chain_output[[3]]$samples[,"a"])))$quantiles[c(1,5)]
b_CI <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,"b"]),
                          mcmc(chain_output[[2]]$samples[,"b"]),
                          mcmc(chain_output[[3]]$samples[,"b"])))$quantiles[c(1,5)]
alpha_CI <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^alpha",col)]),
                              mcmc(chain_output[[2]]$samples[,grep("^alpha",col)]),
                              mcmc(chain_output[[3]]$samples[,grep("^alpha",col)])))$quantiles[,c(1,5)]
eta_CI <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^eta",col)]),
                            mcmc(chain_output[[2]]$samples[,grep("^eta",col)]),
                            mcmc(chain_output[[3]]$samples[,grep("^eta",col)])))$quantiles[,c(1,5)]


## calculate Rhat statistics 
beta1_Rhat <- coda::gelman.diag(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^beta1_ori",col)]),
                                          mcmc(chain_output[[2]]$samples[,grep("^beta1_ori",col)]),
                                          mcmc(chain_output[[3]]$samples[,grep("^beta1_ori",col)])))$psrf[,1]
Psi1_Rhat <- coda::gelman.diag(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^Psi1",col)]),
                                         mcmc(chain_output[[2]]$samples[,grep("^Psi1",col)]),
                                         mcmc(chain_output[[3]]$samples[,grep("^Psi1",col)])))$psrf[,1]
gamma1_Rhat <- coda::gelman.diag(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^gamma1",col)]),
                                           mcmc(chain_output[[2]]$samples[,grep("^gamma1",col)]),
                                           mcmc(chain_output[[3]]$samples[,grep("^gamma1",col)])))$psrf[,1]
sigma_b1_Rhat <- coda::gelman.diag(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^sigma_b1_ori",col)]),
                                             mcmc(chain_output[[2]]$samples[,grep("^sigma_b1_ori",col)]),
                                             mcmc(chain_output[[3]]$samples[,grep("^sigma_b1_ori",col)])))$psrf[,1]
Omega_b1_Rhat <- coda::gelman.diag(mcmc.list(mcmc(Omega_b1_sample[[1]]),
                                             mcmc(Omega_b1_sample[[2]]),
                                             mcmc(Omega_b1_sample[[3]])))$psrf[,1]

beta2_Rhat <- coda::gelman.diag(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^beta2_ori",col)]),
                                          mcmc(chain_output[[2]]$samples[,grep("^beta2_ori",col)]),
                                          mcmc(chain_output[[3]]$samples[,grep("^beta2_ori",col)])))$psrf[,1]
Psi2_Rhat <- coda::gelman.diag(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^Psi2",col)]),
                                         mcmc(chain_output[[2]]$samples[,grep("^Psi2",col)]),
                                         mcmc(chain_output[[3]]$samples[,grep("^Psi2",col)])))$psrf[,1]
gamma2_Rhat <- coda::gelman.diag(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^gamma2",col)]),
                                           mcmc(chain_output[[2]]$samples[,grep("^gamma2",col)]),
                                           mcmc(chain_output[[3]]$samples[,grep("^gamma2",col)])))$psrf[,1]
sigma_b2_Rhat <- coda::gelman.diag(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^sigma_b2_ori",col)]),
                                             mcmc(chain_output[[2]]$samples[,grep("^sigma_b2_ori",col)]),
                                             mcmc(chain_output[[3]]$samples[,grep("^sigma_b2_ori",col)])))$psrf[,1]
Omega_b2_Rhat <- coda::gelman.diag(mcmc.list(mcmc(Omega_b2_sample[[1]]),
                                             mcmc(Omega_b2_sample[[2]]),
                                             mcmc(Omega_b2_sample[[3]])))$psrf[,1]
a_Rhat <- coda::gelman.diag(mcmc.list(mcmc(chain_output[[1]]$samples[,"a"]),
                                      mcmc(chain_output[[2]]$samples[,"a"]),
                                      mcmc(chain_output[[3]]$samples[,"a"])))$psrf[,1]
b_Rhat <- coda::gelman.diag(mcmc.list(mcmc(chain_output[[1]]$samples[,"b"]),
                                      mcmc(chain_output[[2]]$samples[,"b"]),
                                      mcmc(chain_output[[3]]$samples[,"b"])))$psrf[,1]
alpha_Rhat <- coda::gelman.diag(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^alpha",col)]),
                                          mcmc(chain_output[[2]]$samples[,grep("^alpha",col)]),
                                          mcmc(chain_output[[3]]$samples[,grep("^alpha",col)])))$psrf[,1]
eta_Rhat <- coda::gelman.diag(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^eta",col)]),
                                        mcmc(chain_output[[2]]$samples[,grep("^eta",col)]),
                                        mcmc(chain_output[[3]]$samples[,grep("^eta",col)])))$psrf[,1]

summary_table <- data.frame("Estimate" = c(beta1_est, Psi1_est, gamma1_est,
                                           sigma_b1_est, Omega_b1_est,
                                           beta2_est, Psi2_est, gamma2_est,
                                           sigma_b2_est, Omega_b2_est,
                                           a_est, b_est, alpha_est, eta_est),
                            "2.5%" = c(beta1_CI[,1], Psi1_CI[1], gamma1_CI[1],
                                       sigma_b1_CI[,1], Omega_b1_CI[1],
                                       beta2_CI[,1], Psi2_CI[1], gamma2_CI[1],
                                       sigma_b2_CI[,1], Omega_b2_CI[1],
                                       a_CI[1], b_CI[1],
                                       alpha_CI[,1], eta_CI[,1]),
                            "97.5%" = c(beta1_CI[,2], Psi1_CI[2], gamma1_CI[2],
                                        sigma_b1_CI[,2], Omega_b1_CI[2],
                                        beta2_CI[,2], Psi2_CI[2], gamma2_CI[2],
                                        sigma_b2_CI[,2], Omega_b2_CI[2],
                                        a_CI[2], b_CI[2],
                                        alpha_CI[,2], eta_CI[,2]),
                            "Rhat" = c(beta1_Rhat, Psi1_Rhat, gamma1_Rhat,
                                       sigma_b1_Rhat, Omega_b1_Rhat,
                                       beta2_Rhat, Psi2_Rhat, gamma2_Rhat,
                                       sigma_b2_Rhat, Omega_b2_Rhat,
                                       a_Rhat, b_Rhat,
                                       alpha_Rhat, eta_Rhat),
                            row.names = c("beta1[1]", "beta1[2]",
                                          "beta1[3]", "beta1[4]",
                                          "Psi1", "gamma1",
                                          "sigma_b1[1]", "sigma_b1[2]",
                                          "Omega_b1[1,2]",
                                          "beta2[1]", "beta2[2]",
                                          "beta2[3]", 
                                          "Psi2", "gamma2",
                                          "sigma_b2[1]", "sigma_b2[2]",
                                          "Omega_b2[1,2]",
                                          "a", "b",
                                          sapply(1:8, FUN = function(t) paste0("alpha[",t,"]")),
                                          sapply(1:8, FUN = function(t) paste0("eta[",t,"]"))))
print(summary_table, digits = 3)

colnames(summary_table)[2:3] <- c("Lower", "Upper")

save(summary_table, file = "summary_table.Rdata")




