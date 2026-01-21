###### Generation of longitudinal biomarker data and time-to-event outcomes ######

## load packages
library(MASS)
library(ggplot2)
library(data.table)
library(GLMMadaptive)
library(dplyr)
library(statmod)
library(tidyr)
library(ggthemes)

k <- Sys.getenv("SLURM_ARRAY_TASK_ID")
print(k)

set.seed(20251210+as.numeric(k)*100)

## number of individuals
N <- 1000

## set up parameter values for biomarker q=1
### population fixed effects
beta1 <- c(3.0, -0.20, 0.04, -0.001)
### variance component of random effects
sigma_b1 <- c(0.23, 0.05)
Omega_b1 <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
### parameters for residual variability
gamma1 <- -0.95
Psi1 <- 0.45

## set up parameter values for biomarker q=2
### population fixed effects
beta2 <- c(6.6, 0.03, -0.05)
### variance component of random effects
sigma_b2 <- c(0.65, 0.25)
Omega_b2 <- matrix(c(1, 0.18, 0.18, 1), nrow = 2)
### parameters for residual variability
gamma2 <- 0.50
Psi2 <- 1.0

### prior for correlation
a <- 5.3
b <- 12

## parameters for survival model
alpha <- c(3.7, -0.20, 0.15, -0.04,
           0.10, 0.25, 0.10)
eta <- c(-3.8, 0.85, -0.9, -0.02,
         -0.28, -0.36, -0.2)

## generate follow-up time
t <- matrix(0, nrow = N, ncol = 15)
n <- rep(0, N)  ## number of observations per individual
for(i in 1:N){
  drop_i <- min(rpois(1, lambda = 5), 9)  ## drop the minimum between 9 and a Poisson-generated value with a mean of 5
  n[i] <- 15-drop_i
  t[i,1:n[i]] <- sort(runif(n = n[i], min = 1, max = 15))  ## follow-up times are sampled from a uniform distribution and sorted afterwards
}

## generate longitudinal biomarker values
marker1 <- matrix(data = Inf, nrow = N, ncol = 15)
marker2 <- matrix(data = Inf, nrow = N, ncol = 15)
b_rm1 <- matrix(nrow = N, ncol = 2)  ## matrix of random effects for biomarker1
b_rm2 <- matrix(nrow = N, ncol = 2)  ## matrix of random effects for biomarker2
var_sq1 <- rep(0, N)  ## vector of residual variability for biomarker1
var_sq2 <- rep(0, N)  ## vector of residual variability for biomarker2
for(i in 1:N){
  ## simulate random effects from a multivariate Gaussian distribution
  bi1 <- mvrnorm(n = 1, mu = rep(0, 2),
                 Sigma = diag(sigma_b1) %*% Omega_b1 %*% diag(sigma_b1))
  bi2 <- mvrnorm(n = 1, mu = rep(0, 2),
                 Sigma = diag(sigma_b2) %*% Omega_b2 %*% diag(sigma_b2))
  b_rm1[i,] <- bi1
  b_rm2[i,] <- bi2
  
  ## simulate residual variabilities from a log-normal distribution
  var_i1 <- rlnorm(n = 1, meanlog = gamma1, sdlog = Psi1)
  var_i2 <- rlnorm(n = 1, meanlog = gamma2, sdlog = Psi2)
  sdi <- c(sqrt(var_i1), sqrt(var_i2))
  var_sq1[i] <- var_i1
  var_sq2[i] <- var_i2
  
  ## simulate residual correlations from a transformed beta distribution
  ri <- rbeta(n = 1, shape1 = a, shape2 = b)*2-1
  Omega_q <- matrix(c(1, ri, ri, 1), nrow = 2)  ## construct a correlation matrix
  
  ## simulate longitudinal biomarkers from a multivariate-Gaussian distribution
  for(j in 1:n[i]){
    mu_i1 <- beta1[1]+beta1[2]*t[i,j]+beta1[3]*t[i,j]^2+beta1[4]*t[i,j]^3+
      bi1[1]+bi1[2]*t[i,j]
    mu_i2 <- beta2[1]+beta2[2]*t[i,j]+beta2[3]*t[i,j]^2+
      bi2[1]+bi2[2]*t[i,j]
    
    marker_i <- mvrnorm(n = 1, mu = c(mu_i1, mu_i2),
                        Sigma = diag(sdi) %*% Omega_q %*% diag(sdi))
    marker1[i,j] <- marker_i[1]
    marker2[i,j] <- marker_i[2]
  }
}

## set left-censoring rate of marker2 to be 20% (change probs argument to change the proportion of biomarkers below LOD)
LOD <- quantile(marker2[marker2!=Inf], probs = 0.2)
marker_lc <- ifelse(marker2 <= LOD, LOD, marker2)  ## biomarker2 values below LOD were set to LOD

sum(marker_lc==LOD)/sum(n)  ## calculate the proportion of biomarker below LOD
print(LOD)
save(LOD, file = paste0(k, ".result/LOD.Rdata"))

## visualize the simulated longitudinal biomarker data
graph.data <- data.frame()
for(i in 1:N){
  df_i <- data.frame("ID" <- rep(i, n[i]),
                     "times" <- t[i, 1:n[i]],
                     "marker_1" <- marker1[i, 1:n[i]],
                     "marker_2" <- marker2[i, 1:n[i]],
                     "marker2_afterLOD" <- marker_lc[i, 1:n[i]])
  graph.data <- rbind(graph.data, df_i)
}
colnames(graph.data) <- c("ID", "times", "Marker 1", 
                          "Marker 2 before censoring", 
                          "Marker 2 after censoring")

graph.data.long <- graph.data %>%
  pivot_longer(cols = 3:5, names_to = "Biomarker type",
               values_to = "value")

graph.data.long$`Biomarker type` <- factor(graph.data.long$`Biomarker type`,
                                           levels = c("Marker 1", 
                                                      "Marker 2 before censoring", 
                                                      "Marker 2 after censoring"),
                                           ordered = TRUE)

## make a graph of 9 randomly selected individuals
ggplot(data = graph.data.long[graph.data.long$ID %in% sample(1:N, size = 9),])+
  geom_point(aes(x = times, y = value, col = `Biomarker type`, shape = `Biomarker type`))+
  geom_line(aes(x = times, y = value, col = `Biomarker type`))+
  geom_hline(yintercept = -1, linetype = "dashed")+
  facet_wrap(.~ID)+
  scale_shape_manual(values = c(16, 15, 17))+
  scale_color_tableau()+
  theme_bw()+
  ylab("Biomarker values")+
  xlab("Time")+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 9.5),
        legend.key.width = unit(1.5, "cm"))

## simulate time-to-event outcomes with some rate of censoring
## covariates: standardized random effects and residual variances
sim_surv_data <- function(N, b_rm1, b_rm2, var_sq1, var_sq2, alpha, eta,
                         marker1, marker2, marker_lc){
  ## standardize random effects
  ## first biomarker
  b11_standard <- (b_rm1[,1]-0)/sigma_b1[1]
  b12_standard <- (b_rm1[,2]-0)/sigma_b1[2]
  ## second biomarker
  b21_standard <- (b_rm2[,1]-0)/sigma_b2[1]
  b22_standard <- (b_rm2[,2]-0)/sigma_b2[2]
  
  ## standardize residual variability
  ## first biomarker
  var_sq1_standard <- (log(var_sq1)-gamma1)/Psi1
  ## second biomarker
  var_sq2_standard <- (log(var_sq2)-gamma2)/Psi2
  
  ## generate y_0 (initial status)
  y_0 <- rep(0, N)
  for(i in 1:N){
    y_0[i] <- exp(alpha[1]+alpha[2]*b11_standard[i]+
                    alpha[3]*b12_standard[i]+alpha[4]*var_sq1_standard[i]+
                    alpha[5]*b21_standard[i]+
                    alpha[6]*b22_standard[i]+alpha[7]*var_sq2_standard[i])
  }
  ## generate zeta (drift parameter)
  zeta <- rep(0, N)
  for(i in 1:N){
    zeta[i] <- eta[1]+eta[2]*b11_standard[i]+
        eta[3]*b12_standard[i]+eta[4]*var_sq1_standard[i]+
        eta[5]*b21_standard[i]+
        eta[6]*b22_standard[i]+eta[7]*var_sq2_standard[i]
  }

  ## simulate event times from threshold regression (inverse-Gaussian distribution)
  event_time <- rinvgauss(N, mean = -y_0/mu, shape = y_0^2)
  ## simulate censoring times from an independent exponential distribution
  censor_time <- 4+rexp(N, rate = 1/14)
  ## determine observed event times and censoring status
  observe_time <- ifelse(event_time <= censor_time,
                         event_time, censor_time)
  status <- ifelse(event_time <= censor_time,
                   1, 0)
  test_data <- data.frame(observe_time, status, y_0, mu, 
                          b11_standard, b12_standard,
                          var_sq1_standard,
                          b21_standard, b22_standard,
                          var_sq2_standard,
                          marker1, marker2, marker_lc)
  return(test_data)
}

surv_biomarker_data <- sim_surv_data(N = N, b_rm1 = b_rm1, b_rm2 = b_rm2,
                                     var_sq1 = var_sq1, var_sq2 = var_sq2,
                                     alpha = alpha, eta = eta,
                                     marker1 = marker1, marker2 = marker2,
                                     marker_lc = marker_lc)
surv_biomarker_data <- na.omit(surv_data)

save(surv_biomarker_data, file = paste0(k, ".result/surv_biomarker_data.Rdata"))











