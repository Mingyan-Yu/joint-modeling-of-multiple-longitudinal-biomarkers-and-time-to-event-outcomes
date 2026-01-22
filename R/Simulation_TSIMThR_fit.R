###### Two-stage alternative approach: TSIMThR ######

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
load(paste0(k, ".result/surv_data.Rdata"))
load(paste0(k, ".result/t.Rdata")
load(paste0(k, ".result/n.Rdata")
load(paste0(k, ".result/marker1.Rdata")
load(paste0(k, ".result/marker2.Rdata")
load(paste0(k, ".result/marker_lc.Rdata")

N <- nrow(surv_data)

### Stage 1: Fit biomarker-only NIMBLE model to two biomarkers

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

## export the user-defined function to the clusters
clusterExport(cl, "uppertri_mult_diag")
clusterExport(cl, "registerDistributions")
clusterEvalQ(cl, {
  library(nimble)
})

## construct the NIMBLE model for two biomarkers with one subject to detection limits
## and time-to-event outcomes
## with measurement time standardized
Two_biomarker_nimble <- nimbleCode({
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
    b1[i,1:2] ~ dmnorm(mean_b_prior[1:2], cholesky = Ucorr_b1_diag[1:2,1:2], prec_param = 0)
    
    b2[i,1:2] ~ dmnorm(mean_b_prior[1:2], cholesky = Ucorr_b2_diag[1:2,1:2], prec_param = 0)
    
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
  
  ## priors for covariance matrix
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
  }
})

## construct data for fitting biomarker model
## turn left-censored values into NA
marker_lc <- ifelse(marker_lc==LOD, NA, marker_lc)
## construct a censoring indicator for biomarker2
censored_marker <- marker_lc
censored_marker[] <- as.numeric(!is.na(marker_lc))

Two_biomarker_data <- list(marker1 = marker1, marker2 = marker_lc,
                           censored_marker = censored_marker)
save(Two_biomarker_data, file = paste0(k, ".result/Two_biomarker_data.Rdata"))

## standardize measurement time
mean_t = mean(t[t!=0])
sd_t = sd(t[t!=0])
t_std = (t-mean_t)/sd_t

## construct constants for fitting the biomarker model
Two_biomarker_constants <- list(N = N, n = n, time = t_std, LOD = LOD,
                                mean_t = mean_t, sd_t = sd_t)
save(Two_biomarker_constants, file = paste0(k, ".result/Two_biomarker_constants.Rdata"))

## specify some initial values for some parameters
Two_biomarker_inits <- list()

## for beta1, sigma_b1, Omega_b1, gamma1 and Psi1: 
## fit a linear mixed effects model on marker1 by the standardized measurement times
marker1_data <- data.frame()
for(i in 1:N){
  marker1_i <- data.frame("ID" = rep(i, n[i]),
                          "Time" = t_std[i,1:n[i]],
                          "Marker1" = marker1[i,1:n[i]])
  marker1_data <- rbind(marker1_data, marker1_i)
}

marker1_fit <- lme(fixed = Marker1 ~ Time + I(Time^2) + I(Time^3),
                   random = ~ 1 + Time | ID,
                   data = marker1_data)

Two_biomarker_inits$beta1 <- fixef(marker1_fit) ## from fixed effects
Two_biomarker_inits$sigma_b1 <- apply(ranef(marker1_fit), MARGIN = 2, FUN = sd) ## from sd of random effects
Two_biomarker_inits$Ucorr_b1 <- rlkj_corr_cholesky(n = 1, eta = 1, p = 2) ## randomly generate from a lkj distribution
Two_biomarker_inits$Ucov_b1 <- uppertri_mult_diag(Two_biomarker_inits$Ucorr_b1[1:2,1:2], 
                                                  Two_biomarker_inits$sigma_b1[1:2])
marker1_data$residual <- marker1_data$Marker1 - marker1_fit$fitted[,2]  ## obtain individual-specific residuals
## calculate individual-specific residual variance
var_marker1 <- marker1_data %>%
  group_by(ID) %>%
  summarise(Si1 = var(residual))
Two_biomarker_inits$gamma1 <- mean(log(var_marker1$Si1)) ## from mean of log-transformed residual variance
Two_biomarker_inits$Psi1 <- sd(log(var_marker1$Si1)) #from sd of log-transformed residual variance

Two_biomarker_inits$var_sq1 <- rlnorm(Two_biomarker_constants$N, 
                                      meanlog = Two_biomarker_inits$gamma1,
                                      sdlog = Two_biomarker_inits$Psi1)

## for beta2, sigma_b2, Omega_b2, gamma2 and Psi2:
## fit a tobit mixed effects model on marker2 by the standardized measurement times
marker2_data <- data.frame()
for(i in 1:N){
  marker2_i <- data.frame("ID" = rep(i, n[i]),
                          "Time" = t_std[i,1:n[i]],
                          "Marker2" = ifelse(is.na(marker_lc[i,1:n[i]]), LOD,
                                             marker_lc[i,1:n[i]]))
  marker2_data <- rbind(marker2_data, marker2_i)
}
marker2_data$censor <- ifelse(marker2_data$Marker2 == LOD, 1, 0)

marker2_fit <- mixed_model(fixed = cbind(Marker2, censor) ~ Time + I(Time^2),
                           random = ~ 1 + Time | ID,
                           data = marker2_data,
                           family = censored.normal())

Two_biomarker_inits$beta2 <- fixef(marker2_fit) ## from fixed effects
Two_biomarker_inits$sigma_b2 <- apply(ranef(marker2_fit), MARGIN = 2, FUN = sd) ## from sd of random effects
Two_biomarker_inits$Ucorr_b2 <- rlkj_corr_cholesky(n = 1, eta = 1, p = 2) ## randomly generate from a lkj distribution
Two_biomarker_inits$Ucov_b2 <- uppertri_mult_diag(Two_biomarker_inits$Ucorr_b2[1:2,1:2], 
                                                  Two_biomarker_inits$sigma_b2[1:2])
## obtain individual-specific residuals
marker2_data$residual <- marker2_data$Marker2 - GLMMadaptive::fitted(marker2_fit, type = "subject_specific")
## calculate individual-specific residual variance
var_marker2 <- marker2_data %>%
  group_by(ID) %>%
  summarise(Si2 = var(residual))
Two_biomarker_inits$gamma2 <- mean(log(var_marker2$Si2)) ## from mean of log-transformed residual variance
Two_biomarker_inits$Psi2 <- sd(log(var_marker2$Si2)) #from sd of log-transformed residual variance

Two_biomarker_inits$var_sq2 <- rlnorm(Two_biomarker_constants$N, meanlog = Two_biomarker_inits$gamma2,
                                      sdlog = Two_biomarker_inits$Psi2)

## for a and b: calculate the correlation between marker1 residual and marker2 residual
marker_cor_data <- data.frame("ID" = marker1_data$ID,
                              "Marker1_res" = marker1_data$residual,
                              "Marker2_res" = marker2_data$residual)
marker_cor <- marker_cor_data %>%
  group_by(ID) %>%
  summarise(marker_cor = cor(Marker1_res, Marker2_res))
marker_cor_transform <- (marker_cor$marker_cor+1)/2

Two_biomarker_inits$a <- mean(marker_cor_transform, na.rm = TRUE)^2*
  (1-mean(marker_cor_transform, na.rm = TRUE))/var(marker_cor_transform, na.rm = TRUE) - 
  mean(marker_cor_transform, na.rm = TRUE)
Two_biomarker_inits$b <- Two_biomarker_inits$a*(1/mean(marker_cor_transform, na.rm = TRUE)-1)

print(Two_biomarker_inits$a)
print(Two_biomarker_inits$b)

Two_biomarker_inits$t_r <- rbeta(Two_biomarker_constants$N, Two_biomarker_inits$a, Two_biomarker_inits$b)
Two_biomarker_inits$r <- Two_biomarker_inits$t_r*2-1

## specify the rest of the initial values
## for biomarker2 below LOD, specify the initial values for the imputation as LOD-1
Two_biomarker_inits$marker2 <- ifelse(is.na(marker_lc), LOD-1, NA)

save(Two_biomarker_inits, file = paste0(k, ".result/Two_biomarker_inits.Rdata"))

## start to build and fit the NIMBLE model in parallel
run_MCMC_parallel <- function(seed, data, constants, code, inits){
  library(nimble)
  
  ## build the nimble model
  Two_biomarker_model <- nimbleModel(code = code, name = "Two_biomarker_model",
                                     data = data, constants = constants,
                                     inits = inits)
  
  ## configure the model
  Two_biomarker_model.config <- configureMCMC(Two_biomarker_model, print = TRUE)
  Two_biomarker_model.config$addMonitors("biomarker_predictor", "Omega_b1", "Omega_b2",
                                         "beta1_ori", "beta2_ori",
                                         "sigma_b1_ori", "sigma_b2_ori")
  
  ## change samplers for some parameters
  ### replace samplers for beta1 and beta2 with automated factor slice samplers
  Two_biomarker_model.config$removeSamplers("beta1[]", "beta2[]")
  Two_biomarker_model.config$addSampler(target = c("beta1[]", "beta2[]"),
                                        type = "AF_slice")
                                
  ### replace samplers for some other parameters with slice samplers
  Two_biomarker_model.config$removeSamplers(c("sigma_b1[1]", "sigma_b2[1]",
                                              "sigma_b1[2]", "sigma_b2[2]",
                                              "Psi1", "Psi2", "a", "b"))
  Two_biomarker_model.config$addSampler(target = "sigma_b1[]", type = "AF_slice")
  Two_biomarker_model.config$addSampler(target = "sigma_b2[]", type = "AF_slice")
  Two_biomarker_model.config$addSampler(target = "Psi1",
                                        type = "slice")
  Two_biomarker_model.config$addSampler(target = "Psi2",
                                        type = "slice")
  Two_biomarker_model.config$addSampler(target = "a",
                                        type = "slice")
  Two_biomarker_model.config$addSampler(target = "b",
                                        type = "slice")
  
  ## build an MCMC object for this model
  Two_biomarker_model.MCMC <- buildMCMC(Two_biomarker_model.config)
  
  ## compile the model and the MCMC object
  Com_Two_biomarker.Model <- compileNimble(Two_biomarker_model)
  Com_Two_biomarker.MCMC <- compileNimble(Two_biomarker_model.MCMC, project = Two_biomarker_model)
  
  Sys.time()
  ## run the model
  run_Two_biomarker_samples <- runMCMC(Com_Two_biomarker.MCMC, niter = 60000,
                                       nburnin = 30000, progressBar = TRUE,
                                       summary = TRUE, setSeed = seed)
  Sys.time()
  
  return(run_Two_biomarker_samples)
}

chain_output <- parLapply(cl = cl, X = 20260115+as.numeric(k)+(1:3)*100,
                          fun = run_MCMC_parallel,
                          data = Two_biomarker_data,
                          code = Two_biomarker_nimble,
                          constants = Two_biomarker_constants,
                          inits = Two_biomarker_inits)

col <- colnames(chain_output[[1]]$samples)
n.sims <- nrow(chain_output[[1]]$samples)

## post-sampling processing
## calculate the correlation matrix for the random effects
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
pdf(paste0(k, ".result/traceplot_biomarker.pdf"))
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

summary_table <- data.frame("Estimate" = c(beta1_est, Psi1_est, gamma1_est,
                                           sigma_b1_est, Omega_b1_est,
                                           beta2_est, Psi2_est, gamma2_est,
                                           sigma_b2_est, Omega_b2_est,
                                           a_est, b_est),
                            "2.5%" = c(beta1_CI[,1], Psi1_CI[1], gamma1_CI[1],
                                       sigma_b1_CI[,1], Omega_b1_CI[1],
                                       beta2_CI[,1], Psi2_CI[1], gamma2_CI[1],
                                       sigma_b2_CI[,1], Omega_b2_CI[1],
                                       a_CI[1], b_CI[1]),
                            "97.5%" = c(beta1_CI[,2], Psi1_CI[2], gamma1_CI[2],
                                        sigma_b1_CI[,2], Omega_b1_CI[2],
                                        beta2_CI[,2], Psi2_CI[2], gamma2_CI[2],
                                        sigma_b2_CI[,2], Omega_b2_CI[2],
                                        a_CI[2], b_CI[2]),
                            "Rhat" = c(beta1_Rhat, Psi1_Rhat, gamma1_Rhat,
                                       sigma_b1_Rhat, Omega_b1_Rhat,
                                       beta2_Rhat, Psi2_Rhat, gamma2_Rhat,
                                       sigma_b2_Rhat, Omega_b2_Rhat,
                                       a_Rhat, b_Rhat),
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
                                          "a", "b"))
print(summary_table, digits = 3)

colnames(summary_table)[2:3] <- c("Lower", "Upper")

save(summary_table, file = paste0(k, ".result/summary_table_biomarker.Rdata"))

### Stage 2: Fit the threshold regression model on posterior means obtained from the biomarker model fit

## extract matrix of biomarker predictors and calculate posterior means
biomarker_predictors <- array(dim = c(3, N, 7))
for(i in 1:3){
  biomarker_predictors_i <- array(dim = c(N, 7))
  for(j in 1:N){
    for(l in 1:7){
      ## calculate the posterior mean from each chain output
      biomarker_predictors_i[j,l] <- mean(chain_output[[i]]$samples[,paste0("biomarker_predictor[", j, ", ", l, "]")])
    }
  }
  biomarker_predictors[i,,] <- biomarker_predictors_i
}

## calculate the average across posterior means from each chain
biomarker_predictors <- apply(biomarker_predictors, MARGIN = c(2,3), FUN = mean)

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

clusterExport(cl, c("dinversegauss", "rinversegauss"))
clusterExport(cl, "registerDistributions")
clusterEvalQ(cl, library(nimble))

clusterEvalQ(cl, {
  registerDistributions(list(
    dinversegauss = list(
      BUGSdist = "dinversegauss(mean_mu, lambda)",
      types = c("value=double(0)", "mean_mu=double(0)", "lambda=double(0)")
    )
  ))
})

## construct the NIMBLE model for threshold regression
threg_nimble <- nimbleCode({
  for(i in 1:N){
    mean_mu[i] <- -exp(inprod(biomarker_predictor[i,],alpha[1:7]))/inprod(biomarker_predictor[i,],eta[1:7])
    lambda[i] <- exp(inprod(biomarker_predictor[i,],alpha[1:7]))^2
    
    censored[i] ~ dinterval(t[i], c[i])
    t[i] ~ dinversegauss(mean_mu[i], lambda[i])
  }
  
  for(j in 1:7){
    alpha[j] ~ dnorm(0, sd = 10)
    eta[j] ~ dnorm(0, sd = 10)
  }
})

## construct data for fitting threshold regression
threg_data <- list()
threg_data$censored <- as.numeric(surv_data$status==0)
threg_data$t <- ifelse(surv_data$status==1, surv_data$observe_time, NA)

threg_constants <- list()
threg_constants$c <- ifelse(surv_data$status==1, Inf, surv_data$observe_time)
threg_constants$N <- N
threg_constants$biomarker_predictor <- biomarker_predictors

run_MCMC_parallel <- function(seed, data, constants, code){
  library(nimble)
  
  ## initialize parameter values
  alpha_init <- rnorm(n = 7, mean = 0, sd = 0.5)
  eta_init <- rnorm(n = 7, mean = 0, sd = 0.5)
  ## for right-censored event times, specify the initial values for imputation as censoring time + 2
  t_init <- rep(NA, constants$N)
  t_init[data$censored==1] <- constants$c[data$censored==1] + 2
  
  threg_init <- list(alpha = alpha_init, eta = eta_init, t = t_init)
  
  ## build the nimble model
  threg_model <- nimbleModel(code = code, name = "Threshold_regression",
                             data = data, inits = threg_init,
                             constants = constants)
  
  ## configure the model
  threg_model.config <- configureMCMC(threg_model, print = TRUE,
                                      monitors = c("alpha", "eta"))
  
  ## replace samplers for alpha and eta with automated factor slice sampler
  threg_model.config$removeSamplers("alpha[]", "eta[]")
  threg_model.config$addSampler(target = c("alpha[]", "eta[]"),
                                type = "AF_slice")
  
  ## build an MCMC object for this model
  threg_model.MCMC <- buildMCMC(threg_model.config)
  
  ## compile the model and the MCMC object
  Com_threg.Model <- compileNimble(threg_model)
  Com_threg.MCMC <- compileNimble(threg_model.MCMC, project = threg_model)
  
  Sys.time()
  ## run the model
  run_threg_samples <- runMCMC(Com_threg.MCMC, niter = 6000,
                               nburnin = 3000, progressBar = TRUE,
                               summary = TRUE, setSeed = seed)
  Sys.time()
  
  return(run_threg_samples)
}

chain_output <- parLapply(cl = cl, X = 20251208+as.numeric(k)+(1:3)*100,
                          fun = run_MCMC_parallel,
                          data = threg_data,
                          code = threg_nimble,
                          constants = threg_constants)

col <- colnames(chain_output[[1]]$samples)
n.sims <- nrow(chain_output[[1]]$samples)

## post-sampling processing
## make a traceplot of the posterior samples
pdf(paste0(k, ".result/traceplot_survival.pdf"))
coda::traceplot(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^alpha",col)]),
                          mcmc(chain_output[[2]]$samples[,grep("^alpha",col)]),
                          mcmc(chain_output[[3]]$samples[,grep("^alpha",col)])))
coda::traceplot(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^eta",col)]),
                          mcmc(chain_output[[2]]$samples[,grep("^eta",col)]),
                          mcmc(chain_output[[3]]$samples[,grep("^eta",col)])))
dev.off()

## calculate posterior mean and 95% CIs
alpha_est <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^alpha",col)]),
                               mcmc(chain_output[[2]]$samples[,grep("^alpha",col)]),
                               mcmc(chain_output[[3]]$samples[,grep("^alpha",col)])))$statistics[,1]
eta_est <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^eta",col)]),
                             mcmc(chain_output[[2]]$samples[,grep("^eta",col)]),
                             mcmc(chain_output[[3]]$samples[,grep("^eta",col)])))$statistics[,1]

alpha_CI <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^alpha",col)]),
                              mcmc(chain_output[[2]]$samples[,grep("^alpha",col)]),
                              mcmc(chain_output[[3]]$samples[,grep("^alpha",col)])))$quantiles[,c(1,5)]
eta_CI <- summary(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^eta",col)]),
                            mcmc(chain_output[[2]]$samples[,grep("^eta",col)]),
                            mcmc(chain_output[[3]]$samples[,grep("^eta",col)])))$quantiles[,c(1,5)]

## calculate Rhat statistics
alpha_Rhat <- coda::gelman.diag(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^alpha",col)]),
                                          mcmc(chain_output[[2]]$samples[,grep("^alpha",col)]),
                                          mcmc(chain_output[[3]]$samples[,grep("^alpha",col)])))$psrf[,1]
eta_Rhat <- coda::gelman.diag(mcmc.list(mcmc(chain_output[[1]]$samples[,grep("^eta",col)]),
                                        mcmc(chain_output[[2]]$samples[,grep("^eta",col)]),
                                        mcmc(chain_output[[3]]$samples[,grep("^eta",col)])))$psrf[,1]

summary_table <- data.frame("Estimate" = c(alpha_est, eta_est),
                            "2.5%" = c(alpha_CI[,1], eta_CI[,1]),
                            "97.5%" = c(alpha_CI[,2], eta_CI[,2]),
                            "Rhat" = c(alpha_Rhat, eta_Rhat),
                            row.names = c(sapply(1:7, FUN = function(t) paste0("alpha[",t,"]")),
                                          sapply(1:7, FUN = function(t) paste0("eta[",t,"]"))))
print(summary_table, digits = 3)

colnames(summary_table)[2:3] <- c("Lower", "Upper")

save(summary_table, file = paste0(k, ".result/summary_table_survival.Rdata"))






























































