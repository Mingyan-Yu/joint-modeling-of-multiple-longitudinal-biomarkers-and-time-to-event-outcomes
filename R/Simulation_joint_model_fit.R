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

