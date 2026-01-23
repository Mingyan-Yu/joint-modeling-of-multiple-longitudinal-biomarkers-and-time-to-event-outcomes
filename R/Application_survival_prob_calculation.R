###### Calculate survival probabilities for an average individual when varying a certain parameter ######

## load packages
library(nimble)
library(ggplot2)

## import posterior samples
load("chain_output.Rdata")

## import the application dataset
load("swan_restrict4_both.Rdata")

col <- colnames(chain_output[[1]]$samples)
## extract regression coefficients for the survival component
alpha_init <- rbind(chain_output[[1]]$samples[,grepl("^alpha", col)], 
                    chain_output[[2]]$samples[,grepl("^alpha", col)], 
                    chain_output[[3]]$samples[,grepl("^alpha", col)])
eta_drift <- rbind(chain_output[[1]]$samples[,grepl("^eta", col)], 
                   chain_output[[2]]$samples[,grepl("^eta", col)], 
                   chain_output[[3]]$samples[,grepl("^eta", col)])


## obtain posterior means of certain parameters
load("summary_table.Rdata")
row <- rownames(summary_table)
sigma_b1 <- summary_table[grepl("^sigma_b1\\[", row), "Estimate"]
gamma1 <- summary_table[grepl("^gamma1", row), "Estimate"]
Psi1 <- summary_table[grepl("^Psi1", row), "Estimate"]
sigma_b2 <- summary_table[grepl("^sigma_b2\\[", row), "Estimate"]
gamma2 <- summary_table[grepl("^gamma2", row), "Estimate"]
Psi2 <- summary_table[grepl("^Psi2", row), "Estimate"]

## write up a function to calculate y_0i (initial status) and zeta_i (drift parameter)
initial_drift <- function(b11, b12, b21, b22, S1_sq, S2_sq, BMI, alpha, eta){
  ## standardize predictors from biomarkers
  b11_prime <- (b11-0)/sigma_b1[1]  ## random intercept for biomarker1
  b12_prime <- (b12-0)/sigma_b1[2]  ## random slope for biomarker1
  b21_prime <- (b21-0)/sigma_b2[1]  ## random intercept for biomarker2
  b22_prime <- (b22-0)/sigma_b2[2]  ## random slope for biomarker2
  S1_prime <- (log(S1_sq)-gamma1)/Psi1  ## residual variability for biomarker1
  S2_prime <- (log(S2_sq)-gamma2)/Psi2  ## residual variability for biomarker2
  zeta_i <- eta[1]+eta[2]*b11_prime+eta[3]*b12_prime+
    eta[4]*S1_prime+eta[5]*b21_prime+eta[6]*b22_prime+
    eta[7]*S2_prime+eta[8]*BMI
  y0_i <- exp(alpha[1]+alpha[2]*b11_prime+alpha[3]*b12_prime+
                alpha[4]*S1_prime+alpha[5]*b21_prime+alpha[6]*b22_prime+
                alpha[7]*S2_prime+alpha[8]*BMI)
  return(list(zeta_i, y0_i))
}

## write up a stable function to calculate CDF of Inverse-Gaussian distribution
IG_CDF <- function(zeta_i, y0_i, Ti){
  comp1 <- pnorm(-(Ti*zeta_i+y0_i)/sqrt(Ti))
  comp2 <- exp(-2*y0_i*zeta_i)
  comp3 <- pnorm(-(y0_i-Ti*zeta_i)/sqrt(Ti))
  return(comp1+comp2*comp3)
}

## write up a function to calculate pdf of Inverse-Gaussian distribution
IG_pdf <- function(zeta_i, y0_i, Ti){
  comp1 <- y0_i/sqrt(2*pi*Ti^3)
  comp2 <- exp(-(Ti*zeta_i+y0_i)^2/(2*Ti))
  return(comp1*comp2)
}

## CDF graph

## a sequence of possible ages of FMP
time <- seq(from = 45, to = 68, by = 0.1)
k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

### vary log(S1^2) by one standard deviation around the mean: log residual variability of biomarker1
### fix b11, b12, b21, b22 to mean=0 and fix BMI to mean
### fix log(S2^2) to gamma2
### log(S1^2): gamma1-/+0.5*Psi1 
### from 45 to 68

## log(S1^2): gamma1-0.5*Psi1 (half sd below the mean)

init_drift <- mapply(function(alpha_row, eta_row){
initial_drift(0, 0, 0, 0, 
              exp(gamma1-0.5*Psi1), exp(gamma2),
              mean(swanDat_restrict4$BMI),
              alpha_row, eta_row)}, 
as.data.frame(t(alpha_init)), as.data.frame(t(eta_drift)))
## calculate CDF at a certain age: time[k]
CDF <- mapply(function(zeta_i, y0_i){IG_CDF(zeta_i, y0_i, time[k]-42)}, 
              init_drift[1,], init_drift[2,])
## obtain 95% CI of CDF
below_S1_CDF_lower95 <- quantile(CDF, probs = 0.025, na.rm = T)
below_S1_CDF_upper95 <- quantile(CDF, probs = 0.975, na.rm = T)
## obtain posterior mean of CDF
below_S1_CDF <- mean(CDF, na.rm = T)

save(below_S1_CDF_lower95, 
     file = paste0("CDF_S1/",k,"below_S1_CDF_lower95.Rdata"))
save(below_S1_CDF_upper95, 
     file = paste0("CDF_S1/",k,"below_S1_CDF_upper95.Rdata"))
save(below_S1_CDF, 
     file = paste0("CDF_S1/",k,"below_S1_CDF.Rdata"))
save(CDF, file = paste0("CDF_S1/",k,"CDF_below.Rdata"))

## log(S1^2): gamma1+0.5*Psi1 (half sd above the mean)

init_drift <- mapply(function(alpha_row, eta_row){
initial_drift(0, 0, 0, 0,
              exp(gamma1+0.5*Psi1), exp(gamma2),
              mean(swanDat_restrict4$BMI),
              alpha_row, eta_row)}, 
as.data.frame(t(alpha_init)), as.data.frame(t(eta_drift)))
## calculate CDF at a certain age: time[k]
CDF <- mapply(function(zeta_i, y0_i){IG_CDF(zeta_i, y0_i, time[k]-42)}, 
              init_drift[1,], init_drift[2,])
## obtain 95% CI of CDF
above_S1_CDF_lower95 <- quantile(CDF, probs = 0.025, na.rm = T)
above_S1_CDF_upper95 <- quantile(CDF, probs = 0.975, na.rm = T)
## obtain posterior mean of CDF
above_S1_CDF <- mean(CDF, na.rm = T)

save(above_S1_CDF_lower95, 
     file = paste0("CDF_S1/",k,"above_S1_CDF_lower95.Rdata"))
save(above_S1_CDF_upper95, 
     file = paste0("CDF_S1/",k,"above_S1_CDF_upper95.Rdata"))
save(above_S1_CDF, 
     file = paste0("CDF_S1/",k,"above_S1_CDF.Rdata"))
save(CDF, file = paste0("CDF_S1/",k,"CDF_above.Rdata"))

### vary b11 by one standard deviation around the mean: random intercept of biomarker1
### fix b12, b21, b22 to mean=0 and fix BMI to mean
### fix log(S1^2) to gamma1, log(S2^2) to gamma2
### b11: -/+0.5*sigma_b1[1]
### from 45 to 68

## b11: -0.5*sigma_b1[1] (half sd below the mean)

init_drift <- mapply(function(alpha_row, eta_row){
initial_drift(-0.5*sigma_b1[1], 0, 0, 0,
              exp(gamma1), exp(gamma2),
              mean(swanDat_restrict4$BMI),
              alpha_row, eta_row)}, 
as.data.frame(t(alpha_init)), as.data.frame(t(eta_drift)))
## calculate CDF at a certain age: time[k]
CDF <- mapply(function(zeta_i, y0_i){IG_CDF(zeta_i, y0_i, time[k]-42)}, 
              init_drift[1,], init_drift[2,])
## obtain 95% CI of CDF
below_b11_CDF_lower95 <- quantile(CDF, probs = 0.025, na.rm = T)
below_b11_CDF_upper95 <- quantile(CDF, probs = 0.975, na.rm = T)
## obtain posterior mean of CDF
below_b11_CDF <- mean(CDF, na.rm = T)

save(below_b11_CDF_lower95, 
     file = paste0("CDF_b11/",k,"below_b11_CDF_lower95.Rdata"))
save(below_b11_CDF_upper95, 
     file = paste0("CDF_b11/",k,"below_b11_CDF_upper95.Rdata"))
save(below_b11_CDF, 
     file = paste0("CDF_b11/",k,"below_b11_CDF.Rdata"))
save(CDF, file = paste0("CDF_b11/",k,"CDF_below.Rdata"))

## b11: 0.5*sigma_b1[1] (half sd above the mean)

init_drift <- mapply(function(alpha_row, eta_row){
initial_drift(0.5*sigma_b1[1], 0, 0, 0, 
              exp(gamma1), exp(gamma2),
              mean(swanDat_restrict4$BMI),
              alpha_row, eta_row)}, 
as.data.frame(t(alpha_init)), as.data.frame(t(eta_drift)))
CDF <- mapply(function(zeta_i, y0_i){IG_CDF(zeta_i, y0_i, time[k]-42)}, 
              init_drift[1,], init_drift[2,])
above_b11_CDF_lower95 <- quantile(CDF, probs = 0.025, na.rm = T)
above_b11_CDF_upper95 <- quantile(CDF, probs = 0.975, na.rm = T)
above_b11_CDF <- mean(CDF, na.rm = T)

save(above_b11_CDF_lower95, 
     file = paste0("CDF_b11/",k,"above_b11_CDF_lower95.Rdata"))
save(above_b11_CDF_upper95, 
     file = paste0("CDF_b11/",k,"above_b11_CDF_upper95.Rdata"))
save(above_b11_CDF, 
     file = paste0("CDF_b11/",k,"above_b11_CDF.Rdata"))
save(CDF, file = paste0("CDF_b11/",k,"CDF_above.Rdata"))

### vary b12 by one standard deviation around the mean: random slope of biomarker1
### fix b11, b21, b22 to mean=0 and fix BMI to mean
### fix log(S1^2) to gamma1, log(S2^2) to gamma2
### b12: -/+0.5*sigma_b1[2]
### from 45 to 68

## b12: -0.5*sigma_b1[2] (half sd below the mean)

init_drift <- mapply(function(alpha_row, eta_row){
initial_drift(0, -0.5*sigma_b1[2], 0, 0,
              exp(gamma1), exp(gamma2),
              mean(swanDat_restrict4$BMI),
              alpha_row, eta_row)}, 
as.data.frame(t(alpha_init)), as.data.frame(t(eta_drift)))
## calculate CDF at a certain age: time[k]
CDF <- mapply(function(zeta_i, y0_i){IG_CDF(zeta_i, y0_i, time[k]-42)}, 
              init_drift[1,], init_drift[2,])
## obtain 95% CI of CDF
below_b12_CDF_lower95 <- quantile(CDF, probs = 0.025, na.rm = T)
below_b12_CDF_upper95 <- quantile(CDF, probs = 0.975, na.rm = T)
## obtain posterior mean of CDF
below_b12_CDF <- mean(CDF, na.rm = T)

save(below_b12_CDF_lower95, 
     file = paste0("CDF_b12/",k,"below_b12_CDF_lower95.Rdata"))
save(below_b12_CDF_upper95, 
     file = paste0("CDF_b12/",k,"below_b12_CDF_upper95.Rdata"))
save(below_b12_CDF, 
     file = paste0("CDF_b12/",k,"below_b12_CDF.Rdata"))
save(CDF, file = paste0("CDF_b12/",k,"CDF_below.Rdata"))

## b12: 0.5*sigma_b1[2] (half sd above the mean)

init_drift <- mapply(function(alpha_row, eta_row){
initial_drift(0, 0.5*sigma_b1[2], 0, 0,
              exp(gamma1), exp(gamma2),
              mean(swanDat_restrict4$BMI),
              alpha_row, eta_row)}, 
as.data.frame(t(alpha_init)), as.data.frame(t(eta_drift)))
## calculate CDF at a certain age: time[k]
CDF <- mapply(function(zeta_i, y0_i){IG_CDF(zeta_i, y0_i, time[k]-42)}, 
              init_drift[1,], init_drift[2,])
## obtain 95% CI of CDF
above_b12_CDF_lower95 <- quantile(CDF, probs = 0.025, na.rm = T)
above_b12_CDF_upper95 <- quantile(CDF, probs = 0.975, na.rm = T)
## obtain posterior mean of CDF
above_b12_CDF <- mean(CDF, na.rm = T)

save(above_b12_CDF_lower95, 
     file = paste0("CDF_b12/",k,"above_b12_CDF_lower95.Rdata"))
save(above_b12_CDF_upper95, 
     file = paste0("CDF_b12/",k,"above_b12_CDF_upper95.Rdata"))
save(above_b12_CDF, 
     file = paste0("CDF_b12/",k,"above_b12_CDF.Rdata"))
save(CDF, file = paste0("CDF_b12/",k,"CDF_above.Rdata"))


### vary log(S2^2) by one standard deviation around the mean: log residual variability of biomarker2
### fix b11, b12, b21, b22 to mean=0 and fix BMI to mean
### fix log(S1^2) to gamma1
### log(S2^2): gamma2-/+0.5*Psi2
### from 45 to 68

## log(S2^2): gamma2-0.5*Psi2 (half sd below the mean)

init_drift <- mapply(function(alpha_row, eta_row){
initial_drift(0, 0, 0, 0, 
              exp(gamma1), exp(gamma2-0.5*Psi2),
              mean(swanDat_restrict4$BMI),
              alpha_row, eta_row)}, 
as.data.frame(t(alpha_init)), as.data.frame(t(eta_drift)))
## calculate CDF at a certain age: time[k]
CDF <- mapply(function(zeta_i, y0_i){IG_CDF(zeta_i, y0_i, time[k]-42)}, 
              init_drift[1,], init_drift[2,])
## obtain 95% CI of CDF
below_S2_CDF_lower95 <- quantile(CDF, probs = 0.025, na.rm = T)
below_S2_CDF_upper95 <- quantile(CDF, probs = 0.975, na.rm = T)
## obtain posterior mean of CDF
below_S2_CDF <- mean(CDF, na.rm = T)

save(below_S2_CDF_lower95, 
     file = paste0("CDF_S2/",k,"below_S2_CDF_lower95.Rdata"))
save(below_S2_CDF_upper95, 
     file = paste0("CDF_S2/",k,"below_S2_CDF_upper95.Rdata"))
save(below_S2_CDF, 
     file = paste0("CDF_S2/",k,"below_S2_CDF.Rdata"))
save(CDF, file = paste0("CDF_S2/",k,"CDF_below.Rdata"))

## log(S2^2): gamma2+0.5*Psi2 (half sd above the mean)

init_drift <- mapply(function(alpha_row, eta_row){
initial_drift(0, 0, 0, 0,
              exp(gamma1), exp(gamma2+0.5*Psi2),
              mean(swanDat_restrict4$BMI),
              alpha_row, eta_row)}, 
as.data.frame(t(alpha_init)), as.data.frame(t(eta_drift)))
## calculate CDF at a certain age: time[k]
CDF <- mapply(function(zeta_i, y0_i){IG_CDF(zeta_i, y0_i, time[k]-42)}, 
              init_drift[1,], init_drift[2,])
## obtain 95% CI of CDF
above_S2_CDF_lower95 <- quantile(CDF, probs = 0.025, na.rm = T)
above_S2_CDF_upper95 <- quantile(CDF, probs = 0.975, na.rm = T)
## obtain posterior mean of CDF
above_S2_CDF <- mean(CDF, na.rm = T)

save(above_S2_CDF_lower95, 
     file = paste0("CDF_S2/",k,"above_S2_CDF_lower95.Rdata"))
save(above_S2_CDF_upper95, 
     file = paste0("CDF_S2/",k,"above_S2_CDF_upper95.Rdata"))
save(above_S2_CDF, 
     file = paste0("CDF_S2/",k,"above_S2_CDF.Rdata"))
save(CDF, file = paste0("CDF_S2/",k,"CDF_above.Rdata"))

### vary b21 by one standard deviation around the mean: random intercept of biomarker2
### fix b11, b12, b22 to mean=0 and fix BMI to mean
### fix log(S1^2) to gamma1, log(S2^2) to gamma2
### b21: -/+0.5*sigma_b2[1]
### from 45 to 68

## b21: -0.5*sigma_b2[1] (half sd below the mean)

init_drift <- mapply(function(alpha_row, eta_row){
initial_drift(0, 0, -0.5*sigma_b2[1], 0,
              exp(gamma1), exp(gamma2),
              mean(swanDat_restrict4$BMI),
              alpha_row, eta_row)}, 
as.data.frame(t(alpha_init)), as.data.frame(t(eta_drift)))
## calculate CDF at a certain age: time[k]
CDF <- mapply(function(zeta_i, y0_i){IG_CDF(zeta_i, y0_i, time[k]-42)}, 
              init_drift[1,], init_drift[2,])
## obtain 95% CI of CDF
below_b21_CDF_lower95 <- quantile(CDF, probs = 0.025, na.rm = T)
below_b21_CDF_upper95 <- quantile(CDF, probs = 0.975, na.rm = T)
## obtain posterior mean of CDF
below_b21_CDF <- mean(CDF, na.rm = T)

save(below_b21_CDF_lower95, 
     file = paste0("CDF_b21/",k,"below_b21_CDF_lower95.Rdata"))
save(below_b21_CDF_upper95, 
     file = paste0("CDF_b21/",k,"below_b21_CDF_upper95.Rdata"))
save(below_b21_CDF, 
     file = paste0("CDF_b21/",k,"below_b21_CDF.Rdata"))
save(CDF, file = paste0("CDF_b21/",k,"CDF_below.Rdata"))

## b21: 0.5*sigma_b2[1] (half sd above the mean)

init_drift <- mapply(function(alpha_row, eta_row){
initial_drift(0, 0, 0.5*sigma_b2[1], 0, 
              exp(gamma1), exp(gamma2),
              mean(swanDat_restrict4$BMI),
              alpha_row, eta_row)}, 
as.data.frame(t(alpha_init)), as.data.frame(t(eta_drift)))
CDF <- mapply(function(zeta_i, y0_i){IG_CDF(zeta_i, y0_i, time[k]-42)}, 
              init_drift[1,], init_drift[2,])
above_b21_CDF_lower95 <- quantile(CDF, probs = 0.025, na.rm = T)
above_b21_CDF_upper95 <- quantile(CDF, probs = 0.975, na.rm = T)
above_b21_CDF <- mean(CDF, na.rm = T)

save(above_b21_CDF_lower95, 
     file = paste0("CDF_b21/",k,"above_b21_CDF_lower95.Rdata"))
save(above_b21_CDF_upper95, 
     file = paste0("CDF_b21/",k,"above_b21_CDF_upper95.Rdata"))
save(above_b21_CDF, 
     file = paste0("CDF_b21/",k,"above_b21_CDF.Rdata"))
save(CDF, file = paste0("CDF_b21/",k,"CDF_above.Rdata"))

### vary b22 by one standard deviation around the mean: random slope of biomarker2
### fix b11, b12, b21 to mean=0 and fix BMI to mean
### fix log(S1^2) to gamma1, log(S2^2) to gamma2
### b22: -/+0.5*sigma_b2[2]
### from 45 to 68

## b22: -0.5*sigma_b2[2] (half sd below the mean)

init_drift <- mapply(function(alpha_row, eta_row){
initial_drift(0, 0, 0, -0.5*sigma_b2[2],
              exp(gamma1), exp(gamma2),
              mean(swanDat_restrict4$BMI),
              alpha_row, eta_row)}, 
as.data.frame(t(alpha_init)), as.data.frame(t(eta_drift)))
## calculate CDF at a certain age: time[k]
CDF <- mapply(function(zeta_i, y0_i){IG_CDF(zeta_i, y0_i, time[k]-42)}, 
              init_drift[1,], init_drift[2,])
## obtain 95% CI of CDF
below_b22_CDF_lower95 <- quantile(CDF, probs = 0.025, na.rm = T)
below_b22_CDF_upper95 <- quantile(CDF, probs = 0.975, na.rm = T)
## obtain posterior mean of CDF
below_b22_CDF <- mean(CDF, na.rm = T)

save(below_b22_CDF_lower95, 
     file = paste0("CDF_b22/",k,"below_b22_CDF_lower95.Rdata"))
save(below_b22_CDF_upper95, 
     file = paste0("CDF_b22/",k,"below_b22_CDF_upper95.Rdata"))
save(below_b22_CDF, 
     file = paste0("CDF_b22/",k,"below_b22_CDF.Rdata"))
save(CDF, file = paste0("CDF_b22/",k,"CDF_below.Rdata"))

## b22: 0.5*sigma_b2[2] (half sd above the mean)

init_drift <- mapply(function(alpha_row, eta_row){
initial_drift(0, 0, 0, 0.5*sigma_b2[2],
              exp(gamma1), exp(gamma2),
              mean(swanDat_restrict4$BMI),
              alpha_row, eta_row)}, 
as.data.frame(t(alpha_init)), as.data.frame(t(eta_drift)))
## calculate CDF at a certain age: time[k]
CDF <- mapply(function(zeta_i, y0_i){IG_CDF(zeta_i, y0_i, time[k]-42)}, 
              init_drift[1,], init_drift[2,])
## obtain 95% CI of CDF
above_b22_CDF_lower95 <- quantile(CDF, probs = 0.025, na.rm = T)
above_b22_CDF_upper95 <- quantile(CDF, probs = 0.975, na.rm = T)
## obtain posterior mean of CDF
above_b22_CDF <- mean(CDF, na.rm = T)

save(above_b22_CDF_lower95, 
     file = paste0("CDF_b22/",k,"above_b22_CDF_lower95.Rdata"))
save(above_b22_CDF_upper95, 
     file = paste0("CDF_b22/",k,"above_b22_CDF_upper95.Rdata"))
save(above_b12_CDF, 
     file = paste0("CDF_b22/",k,"above_b22_CDF.Rdata"))
save(CDF, file = paste0("CDF_b22/",k,"CDF_above.Rdata"))

### vary BMI by one standard deviation around the mean
### fix b11, b12, b21 and b22 to mean=0
### fix log(S1^2) to gamma1, log(S2^2) to gamma2
### BMI: mean-/+2.5
### from 45 to 68

## BMI: mean-0.5sd

init_drift <- mapply(function(alpha_row, eta_row){
  initial_drift(0, 0, 0, 0,
                exp(gamma1), exp(gamma2),
                mean(swanDat_restrict4$BMI)-0.5*sd(swanDat_restrict4$BMI),
                alpha_row, eta_row)}, 
  as.data.frame(t(alpha_init)), as.data.frame(t(eta_drift)))
## calculate CDF at a certain age: time[k]
CDF <- mapply(function(zeta_i, y0_i){IG_CDF(zeta_i, y0_i, time[k]-42)}, 
              init_drift[1,], init_drift[2,])
## obtain 95% CI of CDF
below_BMI_CDF_lower95 <- quantile(CDF, probs = 0.025, na.rm = T)
below_BMI_CDF_upper95 <- quantile(CDF, probs = 0.975, na.rm = T)
## obtain posterior mean of CDF
below_BMI_CDF <- mean(CDF, na.rm = T)

save(below_BMI_CDF_lower95, 
     file = paste0("CDF_BMI/",k,"below_BMI_CDF_lower95.Rdata"))
save(below_BMI_CDF_upper95, 
     file = paste0("CDF_BMI/",k,"below_BMI_CDF_upper95.Rdata"))
save(below_BMI_CDF, 
     file = paste0("CDF_BMI/",k,"below_BMI_CDF.Rdata"))
save(CDF, file = paste0("CDF_BMI/",k,"CDF_below.Rdata"))

## BMI: mean+0.5sd

init_drift <- mapply(function(alpha_row, eta_row){
  initial_drift(0, 0, 0, 0,
                exp(gamma1), exp(gamma2),
                mean(swanDat_restrict4$BMI)+0.5*sd(swanDat_restrict4$BMI),
                alpha_row, eta_row)}, 
  as.data.frame(t(alpha_init)), as.data.frame(t(eta_drift)))
## calculate CDF at a certain age: time[k]
CDF <- mapply(function(zeta_i, y0_i){IG_CDF(zeta_i, y0_i, time[k]-42)}, 
              init_drift[1,], init_drift[2,])
## obtain 95% CI of CDF
above_BMI_CDF_lower95 <- quantile(CDF, probs = 0.025, na.rm = T)
above_BMI_CDF_upper95 <- quantile(CDF, probs = 0.975, na.rm = T)
## obtain posterior mean of CDF
above_BMI_CDF <- mean(CDF, na.rm = T)

save(above_BMI_CDF_lower95, 
     file = paste0("CDF_BMI/",k,"above_BMI_CDF_lower95.Rdata"))
save(above_BMI_CDF_upper95, 
     file = paste0("CDF_BMI/",k,"above_BMI_CDF_upper95.Rdata"))
save(above_BMI_CDF, 
     file = paste0("CDF_BMI/",k,"above_BMI_CDF.Rdata"))
save(CDF, file = paste0("CDF_BMI/",k,"CDF_above.Rdata"))




























