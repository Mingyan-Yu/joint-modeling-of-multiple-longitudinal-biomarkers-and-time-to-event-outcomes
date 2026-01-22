###### Two-stage alternative approach: TSMEThR ######

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

### Stage 1: Fit mixed effects models to two biomarkers
## fit a linear mixed effects model on biomarker1
data_marker1 <- data.frame()
for(i in 1:N){
  data_marker1_i <- data.frame("ID" = i,
                               "time" = t[i, 1:n[i]],
                               "y1" = marker1[i, 1:n[i]])
  data_marker1 <- rbind(data_marker1, data_marker1_i)
}

linear.fit <- lme(fixed = y1 ~ time + I(time^2) + I(time^3),
                  random = ~ 1 + time | ID,
                  data = data_marker1)
summary(linear.fit)

## obtain individual-specific fitted values
data_marker1$fitted <- linear.fit$fitted[,2]
## calculate individual-specific residual variances
Si1 <- data_marker1 %>%
  group_by(ID) %>%
  summarise("Si1" = var(y1-fitted))

## fit a tobit mixed effects model on biomarker2
data_lc <- data.frame()
for(i in 1:N){
  data_lc_i <- data.frame("ID" = i,
                          "time" = t[i,1:n[i]],
                          "y2" = marker_lc[i,1:n[i]])
  data_lc <- rbind(data_lc, data_lc_i)
}
data_lc$censor <- ifelse(data_lc$y2 == LOD, 1, 0)

tobit.fit <- mixed_model(fixed = cbind(y2, censor) ~ time + I(time^2),
                         random = ~ 1 + time | ID,
                         data = data_lc,
                         family = censored.normal())

summary(tobit.fit)

## obtain individual-specific fitted values
data_lc$fitted <- GLMMadaptive::fitted(tobit.fit, type = "subject_specific")
## calculate individual-specific residual variances
Si2 <- data_lc %>%
  group_by(ID) %>%
  summarise("Si2" = var(fitted - y2))

### Stage 2: Fit the threshold regression model on predictors from the longitudinal biomarkers

biomarker_predictors <- matrix(data = NA, nrow = N, ncol = 7)

biomarker_predictors[,1] <- 1
biomarker_predictors[,2] <- (linear.fit$coefficients$random$ID[,1]-mean(linear.fit$coefficients$random$ID[,1]))/sd(linear.fit$coefficients$random$ID[,1])
biomarker_predictors[,3] <- (linear.fit$coefficients$random$ID[,2]-mean(linear.fit$coefficients$random$ID[,2]))/sd(linear.fit$coefficients$random$ID[,2])
biomarker_predictors[,4] <- (log(Si1$Si1)-mean(log(Si1$Si1)))/sd(log(Si1$Si1))
biomarker_predictors[,5] <- (GLMMadaptive::ranef(tobit.fit)[,1]-mean(GLMMadaptive::ranef(tobit.fit)[,1]))/sd(GLMMadaptive::ranef(tobit.fit)[,1])
biomarker_predictors[,6] <- (GLMMadaptive::ranef(tobit.fit)[,2]-mean(GLMMadaptive::ranef(tobit.fit)[,2]))/sd(GLMMadaptive::ranef(tobit.fit)[,2])
biomarker_predictors[,7] <- (log(Si2$Si2)-mean(log(Si2$Si2)))/sd(log(Si2$Si2))

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
## turn right-censored event times into NA
threg_data$t <- ifelse(surv_data$status==1, surv_data$observe_time, NA)

## construct constants for fitting threshold regression
threg_constants <- list()
threg_constants$c <- ifelse(surv_data$status==1, Inf, surv_data$observe_time)
threg_constants$N <- N
threg_constants$biomarker_predictor <- biomarker_predictors

## start to build and fit the NIMBLE model in parallel
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
pdf(paste0(k, ".result/traceplot.pdf"))
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

save(summary_table, file = paste0(k, ".result/summary_table.Rdata"))
























     
