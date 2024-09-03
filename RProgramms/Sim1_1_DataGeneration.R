################################################################################
######                 Data generation of Simulation 1                    ######
################################################################################

# Setting working directory (where functions and Simulation.seed are saved)
setwd("") 

#---------------------- Needed Packages and Functions --------------------------
library(foreach)
library(doParallel)

source("Sim1_0_Functions_DataGeneration.R")

#-------------------- Same parameter in all scenarios  -------------------------
# Number of used cores for Simulations
num.cl <- 35

# Number of iterations for each sub-scenario
n.sim <- 10000 

# Seed for the simulation:
#set.seed(1234567890)
#seed <- sample(x = 1:1000000000, size = n.sim, replace = FALSE)
#save(seed, file = "Sim1_Simulation.seed.RData")
load("Sim1_Simulation.seed.RData")


# Path for saving generated data:
path_data <- "" # e.g. ".../Simulation1/Data/" 


#------------------------- Standard Scenario 1 ---------------------------------
#### Parameters
path <- paste0(path_data, "Scen1")

med.C        <- c(6, 12, 18, 24, 30) 
designHR     <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
p.C          <- c(0.2, 0.4, 0.6) 
p.C.C        <- NA
p.C.T        <- NA
shape.C.T    <- NA

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, p.C, p.C.C, p.C.T, r, med.C, shape.C.T, HR.var, 
                          designHR, n.sim)
# Adding colnames to parameter matrix
colnames(parameters) <- c("beta", "p.C", "p.C.C", "p.C.T", "r", "med.C", 
                          "shape.C.T", "HR.var", "designHR", "n.sim")

# Adding FU parameters, which are dependent on med.C
parameters$accrual.time <- 24 # 2 years = 24 months 
parameters$fu.time      <- 2*parameters$med.C

# Calculating sample size for each sub-scenario and adding n.sim to the 
# parameter matrix
parameters$n.obs.T <- NA
parameters$n.obs.C <- NA
for (i in 1:dim(parameters)[1]) {
  n <- func.sample.size(
    alpha        = alpha,
    beta         = parameters$beta[i],
    r            = parameters$r[i],
    designHR     = parameters$designHR[i],
    med.C        = parameters$med.C[i],
    fu.time      = parameters$fu.time[i],
    accrual.time = parameters$accrual.time[i],
    p.C          = parameters$p.C[i])
  parameters$n.obs.T[i] = n$n.obs.T
  parameters$n.obs.C[i] = n$n.obs.C
}

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path, "..parameters.rds"))


#### Data Generation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), .combine = rbind, 
                  .packages = c("survival")) %dopar% {
                    # Generating Data "n.sim" times for each sub-scenario
                    data.scen.i <- rep(list(NA), para$n.sim)
                    for (i in 1:para$n.sim) {
                      data.scen.i[[i]] <- 
                        DataEXP(
                          n.obs.T        = para$n.obs.T, 
                          n.obs.C        = para$n.obs.C,
                          median.control = para$med.C,
                          accrual.time   = para$accrual.time,
                          fu.time        = para$fu.time,
                          designHR       = para$designHR, 
                          HR.var         = para$HR.var, 
                          p.C            = para$p.C,
                          p.C.C          = para$p.C.C,
                          p.C.T          = para$p.C.T,
                          index.seed     = seed[i])
                    }
                    
                    # Saving data set
                    saveRDS(data.scen.i,
                            file = paste0(path,
                                          "..beta.",     para$beta,
                                          ".p.C.",       para$p.C,
                                          ".p.C.C.",     para$p.C.C,
                                          ".p.C.T.",     para$p.C.T,
                                          ".r.",         para$r,
                                          ".med.C.",     para$med.C,
                                          ".shape.C.T.", para$shape.C.T,
                                          ".HR.var.",    para$HR.var,
                                          ".designHR.",  para$designHR,
                                          ".accrual.",   para$accrual.time,
                                          ".fu.",        para$fu.time,
                                          ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

(time.1 <- Sys.time() - start.time.1)


#------------------------------ Scenario 2 -------------------------------------
#---- Influence of incorrectly assumed designHR for sample size calculation ----

#### Parameters
path <- paste0(path_data, "Scen2")

med.C        <- c(6, 12, 18, 24, 30) 
designHR     <- seq(0.3, 0.9, by = 0.02)
HR.var       <- c(0.8, 0.9, 1.1, 1.2)
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
p.C          <- c(0.2, 0.4, 0.6) 
p.C.C        <- NA
p.C.T        <- NA
shape.C.T    <- NA 

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, p.C, p.C.C, p.C.T, r, med.C, shape.C.T, HR.var, 
                          designHR, n.sim)

# Adding colnames to parameter matrix
colnames(parameters) <- c("beta", "p.C", "p.C.C", "p.C.T", "r", "med.C", 
                          "shape.C.T", "HR.var", "designHR", "n.sim")

# Adding FU parameters, which are dependent on med.C
parameters$accrual.time <- 24 # 2 years = 24 months 
parameters$fu.time      <- 2*parameters$med.C

# Calculating sample size for each sub-scenario and adding it to the parameter 
# matrix
parameters$n.obs.T <- NA
parameters$n.obs.C <- NA
for (i in 1:dim(parameters)[1]) {
  n <- func.sample.size(
    alpha        = alpha,
    beta         = parameters$beta[i],
    r            = parameters$r[i],
    designHR     = parameters$designHR[i],
    med.C        = parameters$med.C[i],
    fu.time      = parameters$fu.time[i],
    accrual.time = parameters$accrual.time[i],
    p.C          = parameters$p.C[i])
  parameters$n.obs.T[i] = n$n.obs.T
  parameters$n.obs.C[i] = n$n.obs.C
}

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path, "..parameters.rds"))

#### Data Generation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), .combine = rbind,
                  .packages = c("survival")) %dopar% {
                    # Generating Data "n.sim" times for each sub-scenario
                    data.scen.i <- rep(list(NA), para$n.sim)
                    for (i in 1:para$n.sim) {
                      data.scen.i[[i]] <- 
                        DataEXP(
                          n.obs.T        = para$n.obs.T, 
                          n.obs.C        = para$n.obs.C,
                          median.control = para$med.C,
                          accrual.time   = para$accrual.time,
                          fu.time        = para$fu.time,
                          designHR       = para$designHR, 
                          HR.var         = para$HR.var, 
                          p.C            = para$p.C,
                          p.C.C          = para$p.C.C,
                          p.C.T          = para$p.C.T,
                          index.seed     = seed[i])
                    }
                    
                    # Saving data set
                    saveRDS(data.scen.i,
                            file = paste0(path,
                                          "..beta.",     para$beta,
                                          ".p.C.",       para$p.C,
                                          ".p.C.C.",     para$p.C.C,
                                          ".p.C.T.",     para$p.C.T,
                                          ".r.",         para$r,
                                          ".med.C.",     para$med.C,
                                          ".shape.C.T.", para$shape.C.T,
                                          ".HR.var.",    para$HR.var,
                                          ".designHR.",  para$designHR,
                                          ".accrual.",   para$accrual.time,
                                          ".fu.",        para$fu.time,
                                          ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

(time.1 <- Sys.time() - start.time.1)


#------------------------------ Scenario 3 -------------------------------------
#--- Influence of different underlying failure time distributions (WEIBULL) ----
#### Parameters
path <- paste0(path_data, "Scen3WEIB")

med.C        <- c(6, 12, 18, 24, 30) 
designHR     <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
p.C          <- c(0.2, 0.4, 0.6) 
p.C.C        <- NA
p.C.T        <- NA
shape.C.T    <- c(0.5, 1.5) 

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, p.C, p.C.C, p.C.T, r, med.C, shape.C.T, HR.var, 
                          designHR, n.sim)

# Adding colnames to parameter matrix
colnames(parameters) <- c("beta", "p.C", "p.C.C", "p.C.T", "r", "med.C", 
                          "shape.C.T", "HR.var", "designHR", "n.sim")

# Adding FU parameters, which are dependent on med.C
parameters$accrual.time <- 24 # 2 years = 24 months 
parameters$fu.time      <- 2*parameters$med.C

# Calculating sample size for each sub-scenario and adding it to the parameter
# matrix
parameters$n.obs.T <- NA
parameters$n.obs.C <- NA
for (i in 1:dim(parameters)[1]) {
  n <- func.sample.size(
    alpha        = alpha,
    beta         = parameters$beta[i],
    r            = parameters$r[i],
    designHR     = parameters$designHR[i],
    med.C        = parameters$med.C[i],
    fu.time      = parameters$fu.time[i],
    accrual.time = parameters$accrual.time[i],
    p.C          = parameters$p.C[i])
  parameters$n.obs.T[i] = n$n.obs.T
  parameters$n.obs.C[i] = n$n.obs.C
}

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path, "..parameters.rds"))

#### Data Generation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), .combine = rbind, 
                  .packages = c("survival")) %dopar% {
                    # Generating Data "n.sim" times for each sub-scenario
                    data.scen.i <- rep(list(NA), para$n.sim)
                    for (i in 1:para$n.sim) {
                      data.scen.i[[i]] <- 
                        DataWEIB(
                          n.obs.T        = para$n.obs.T, 
                          n.obs.C        = para$n.obs.C,
                          median.control = para$med.C,
                          accrual.time   = para$accrual.time,
                          fu.time        = para$fu.time,
                          designHR       = para$designHR, 
                          HR.var         = para$HR.var, 
                          p.C            = para$p.C,
                          shape.C.T      = para$shape.C.T,
                          index.seed     = seed[i])
                    }
                    
                    # Saving data set
                    saveRDS(data.scen.i,
                            file = paste0(path,
                                          "..beta.",     para$beta,
                                          ".p.C.",       para$p.C,
                                          ".p.C.C.",     para$p.C.C,
                                          ".p.C.T.",     para$p.C.T,
                                          ".r.",         para$r,
                                          ".med.C.",     para$med.C,
                                          ".shape.C.T.", para$shape.C.T,
                                          ".HR.var.",    para$HR.var,
                                          ".designHR.",  para$designHR,
                                          ".accrual.",   para$accrual.time,
                                          ".fu.",        para$fu.time,
                                          ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

(time.1 <- Sys.time() - start.time.1)


#------------------------------ Scenario 3 -------------------------------------
#--- Influence of different underlying failure time distributions (GOMPERTZ) ---
#### Parameters
path <- paste0(path_data, "Scen3GOMP")

med.C        <- c(6, 12, 18, 24, 30) 
designHR     <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
p.C          <- c(0.2, 0.4, 0.6) 
p.C.C        <- NA
p.C.T        <- NA
shape.C.T    <- c(-0.2, 0.2) 

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, p.C, p.C.C, p.C.T, r, med.C, shape.C.T, HR.var, 
                          designHR, n.sim)

# Adding colnames to parameter matrix
colnames(parameters) <- c("beta", "p.C", "p.C.C", "p.C.T", "r", "med.C", 
                          "shape.C.T", "HR.var", "designHR", "n.sim")

# Adding FU parameters, which are dependent on med.C
parameters$accrual.time <- 24 # 2 years = 24 months 
parameters$fu.time      <- 2*parameters$med.C

# Calculating sample size for each sub-scenario and adding it to the parameter
# matrix
parameters$n.obs.T <- NA
parameters$n.obs.C <- NA
for (i in 1:dim(parameters)[1]) {
  n <- func.sample.size(
    alpha        = alpha,
    beta         = parameters$beta[i],
    r            = parameters$r[i],
    designHR     = parameters$designHR[i],
    med.C        = parameters$med.C[i],
    fu.time      = parameters$fu.time[i],
    accrual.time = parameters$accrual.time[i],
    p.C          = parameters$p.C[i])
  parameters$n.obs.T[i] = n$n.obs.T
  parameters$n.obs.C[i] = n$n.obs.C
}

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path, "..parameters.rds"))

#### Data Generation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), .combine = rbind, 
                  .packages = c("survival", "flexsurv")) %dopar% {
                    # Generating Data "n.sim" times for each sub-scenario
                    data.scen.i <- rep(list(NA), para$n.sim)
                    for (i in 1:para$n.sim) {
                      data.scen.i[[i]] <- 
                        DataGOMP(
                          n.obs.T        = para$n.obs.T, 
                          n.obs.C        = para$n.obs.C,
                          median.control = para$med.C,
                          accrual.time   = para$accrual.time,
                          fu.time        = para$fu.time,
                          designHR       = para$designHR, 
                          HR.var         = para$HR.var, 
                          p.C            = para$p.C,
                          shape.C.T      = para$shape.C.T,
                          index.seed     = seed[i])
                    }
                    
                    # Saving data set
                    saveRDS(data.scen.i,
                            file = paste0(path,
                                          "..beta.",     para$beta,
                                          ".p.C.",       para$p.C,
                                          ".p.C.C.",     para$p.C.C,
                                          ".p.C.T.",     para$p.C.T,
                                          ".r.",         para$r,
                                          ".med.C.",     para$med.C,
                                          ".shape.C.T.", para$shape.C.T,
                                          ".HR.var.",    para$HR.var,
                                          ".designHR.",  para$designHR,
                                          ".accrual.",   para$accrual.time,
                                          ".fu.",        para$fu.time,
                                          ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

(time.1 <- Sys.time() - start.time.1)


#------------------------------ Scenario 4 -------------------------------------
#--- Influence of non-proportional hazards using late treatment effects for ----
#-------------------------- the treatment group --------------------------------
#### Parameters
path <- paste0(path_data, "Scen4")

med.C        <- c(6, 12, 18, 24, 30) 
designHR     <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
p.C          <- c(0.2, 0.4, 0.6) 
p.C.C        <- NA
p.C.T        <- NA
shape.C.T    <- NA
# treat.effect.start is set to 1/3*med.C due to ESMOs design of using "gain" (med.T-med.C)

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, p.C, p.C.C, p.C.T, r, med.C, shape.C.T, HR.var, 
                          designHR, n.sim)

# Adding colnames to parameter matrix
colnames(parameters) <- c("beta", "p.C", "p.C.C", "p.C.T", "r", "med.C", 
                          "shape.C.T", "HR.var", "designHR", "n.sim")

# Adding FU parameters, which are dependent on med.C
parameters$accrual.time <- 24 # 2 years = 24 months 
parameters$fu.time      <- 2*parameters$med.C

# Calculating sample size for each sub-scenario and adding it to the parameter 
# matrix
parameters$n.obs.T <- NA
parameters$n.obs.C <- NA
for (i in 1:dim(parameters)[1]) {
  n <- func.sample.size(
    alpha        = alpha,
    beta         = parameters$beta[i],
    r            = parameters$r[i],
    designHR     = parameters$designHR[i],
    med.C        = parameters$med.C[i],
    fu.time      = parameters$fu.time[i],
    accrual.time = parameters$accrual.time[i],
    p.C          = parameters$p.C[i])
  parameters$n.obs.T[i] = n$n.obs.T
  parameters$n.obs.C[i] = n$n.obs.C
}

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path, "..parameters.rds"))


#### Data Generation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), .combine = rbind, .packages = c("survival")) %dopar% {
  # Generating Data "n.sim" times for each sub-scenario
  data.scen.i <- rep(list(NA), para$n.sim)
  for (i in 1:para$n.sim) {
    data.scen.i[[i]] <- 
      DataEXPNonProp(
        n.obs.T        = para$n.obs.T, 
        n.obs.C        = para$n.obs.C,
        median.control = para$med.C,
        accrual.time   = para$accrual.time,
        fu.time        = para$fu.time,
        designHR       = para$designHR, 
        HR.var         = para$HR.var, 
        effect.start.T = 1/3 * para$med.C,
        p.C            = para$p.C,
        index.seed     = seed[i])
  }
  
  # Saving data set
  saveRDS(data.scen.i,
          file = paste0(path,
                        "..beta.",     para$beta,
                        ".p.C.",       para$p.C,
                        ".p.C.C.",     para$p.C.C,
                        ".p.C.T.",     para$p.C.T,
                        ".r.",         para$r,
                        ".med.C.",     para$med.C,
                        ".shape.C.T.", para$shape.C.T,
                        ".HR.var.",    para$HR.var,
                        ".designHR.",  para$designHR,
                        ".accrual.",   para$accrual.time,
                        ".fu.",        para$fu.time,
                        ".rds"))
}
# Stopping parallel computing 
stopCluster(cl)

(time.1 <- Sys.time() - start.time.1)


#------------------------------ Scenario 5 -------------------------------------
#---------------------- Influence of unequal sample sizes ----------------------
#### Parameters
path <- paste0(path_data, "Scen5")

med.C        <- c(6, 12, 18, 24, 30) 
designHR     <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- c(0.5, 2)
p.C          <- c(0.2, 0.4, 0.6) 
p.C.C        <- NA
p.C.T        <- NA
shape.C.T    <- NA

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, p.C, p.C.C, p.C.T, r, med.C, shape.C.T, HR.var, 
                          designHR, n.sim)
# Adding colnames to parameter matrix
colnames(parameters) <- c("beta", "p.C", "p.C.C", "p.C.T", "r", "med.C", 
                          "shape.C.T", "HR.var", "designHR", "n.sim")

# Adding FU parameters, which are dependent on med.C
parameters$accrual.time <- 24 # 2 years = 24 months 
parameters$fu.time      <- 2*parameters$med.C

# Calculating sample size for each sub-scenario and adding n.sim to the 
# parameter matrix
parameters$n.obs.T <- NA
parameters$n.obs.C <- NA
for (i in 1:dim(parameters)[1]) {
  n <- func.sample.size(
    alpha        = alpha,
    beta         = parameters$beta[i],
    r            = parameters$r[i],
    designHR     = parameters$designHR[i],
    med.C        = parameters$med.C[i],
    fu.time      = parameters$fu.time[i],
    accrual.time = parameters$accrual.time[i],
    p.C          = parameters$p.C[i])
  parameters$n.obs.T[i] = n$n.obs.T
  parameters$n.obs.C[i] = n$n.obs.C
}

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path, "..parameters.rds"))


#### Data Generation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), .combine = rbind, 
                  .packages = c("survival")) %dopar% {
                    # Generating Data "n.sim" times for each sub-scenario
                    data.scen.i <- rep(list(NA), para$n.sim)
                    for (i in 1:para$n.sim) {
                      data.scen.i[[i]] <- 
                        DataEXP(
                          n.obs.T        = para$n.obs.T, 
                          n.obs.C        = para$n.obs.C,
                          median.control = para$med.C,
                          accrual.time   = para$accrual.time,
                          fu.time        = para$fu.time,
                          designHR       = para$designHR, 
                          HR.var         = para$HR.var, 
                          p.C            = para$p.C,
                          p.C.C          = para$p.C.C,
                          p.C.T          = para$p.C.T,
                          index.seed     = seed[i])
                    }
                    
                    # Saving data set
                    saveRDS(data.scen.i,
                            file = paste0(path,
                                          "..beta.",     para$beta,
                                          ".p.C.",       para$p.C,
                                          ".p.C.C.",     para$p.C.C,
                                          ".p.C.T.",     para$p.C.T,
                                          ".r.",         para$r,
                                          ".med.C.",     para$med.C,
                                          ".shape.C.T.", para$shape.C.T,
                                          ".HR.var.",    para$HR.var,
                                          ".designHR.",  para$designHR,
                                          ".accrual.",   para$accrual.time,
                                          ".fu.",        para$fu.time,
                                          ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

(time.1 <- Sys.time() - start.time.1)


#------------------------------ Scenario 6 -------------------------------------
#--- Influence of using only exponential distributed censoring distribution ----
#---------------------- without administrative censoring -----------------------
#### Parameters
path <- paste0(path_data, "Scen6")

med.C        <- c(6, 12, 18, 24, 30) 
designHR     <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
p.C          <- c(0.2, 0.4, 0.6) 
p.C.C        <- NA
p.C.T        <- NA
shape.C.T    <- NA

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, p.C, p.C.C, p.C.T, r, med.C, shape.C.T, HR.var, 
                          designHR, n.sim)
# Adding colnames to parameter matrix
colnames(parameters) <- c("beta", "p.C", "p.C.C", "p.C.T", "r", "med.C", 
                          "shape.C.T", "HR.var", "designHR", "n.sim")

# Adding FU parameters, which are dependent on med.C
parameters$accrual.time <- NA
parameters$fu.time      <- NA

# Calculating sample size for each sub-scenario and adding n.sim to the 
# parameter matrix
parameters$n.obs.T <- NA
parameters$n.obs.C <- NA
for (i in 1:dim(parameters)[1]) {
  n <- func.sample.size(
    alpha        = alpha,
    beta         = parameters$beta[i],
    r            = parameters$r[i],
    designHR     = parameters$designHR[i],
    med.C        = parameters$med.C[i],
    fu.time      = parameters$fu.time[i],
    accrual.time = parameters$accrual.time[i],
    p.C          = parameters$p.C[i])
  parameters$n.obs.T[i] = n$n.obs.T
  parameters$n.obs.C[i] = n$n.obs.C
}

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path, "..parameters.rds"))


#### Data Generation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), .combine = rbind, 
                  .packages = c("survival")) %dopar% {
                    # Generating Data "n.sim" times for each sub-scenario
                    data.scen.i <- rep(list(NA), para$n.sim)
                    for (i in 1:para$n.sim) {
                      data.scen.i[[i]] <- 
                        DataEXP(
                          n.obs.T        = para$n.obs.T, 
                          n.obs.C        = para$n.obs.C,
                          median.control = para$med.C,
                          accrual.time   = para$accrual.time,
                          fu.time        = para$fu.time,
                          designHR       = para$designHR, 
                          HR.var         = para$HR.var, 
                          p.C            = para$p.C,
                          p.C.C          = para$p.C.C,
                          p.C.T          = para$p.C.T,
                          index.seed     = seed[i])
                    }
                    
                    # Saving data set
                    saveRDS(data.scen.i,
                            file = paste0(path,
                                          "..beta.",     para$beta,
                                          ".p.C.",       para$p.C,
                                          ".p.C.C.",     para$p.C.C,
                                          ".p.C.T.",     para$p.C.T,
                                          ".r.",         para$r,
                                          ".med.C.",     para$med.C,
                                          ".shape.C.T.", para$shape.C.T,
                                          ".HR.var.",    para$HR.var,
                                          ".designHR.",  para$designHR,
                                          ".accrual.",   para$accrual.time,
                                          ".fu.",        para$fu.time,
                                          ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

(time.1 <- Sys.time() - start.time.1)


#------------------------------ Scenario 7 -------------------------------------
#------------ Influence of informative censoring due to treatment --------------
#### Parameters
path <- paste0(path_data, "Scen7")

med.C        <- c(6, 12, 18, 24, 30) 
designHR     <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
p.C          <- NA
p.C.C        <- c(0.2, 0.4)
p.C.T        <- c(0.2, 0.4)
shape.C.T    <- NA

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, p.C, p.C.C, p.C.T, r, med.C, shape.C.T, HR.var, 
                          designHR, n.sim)
# Adding colnames to parameter matrix
colnames(parameters) <- c("beta", "p.C", "p.C.C", "p.C.T", "r", "med.C", 
                          "shape.C.T", "HR.var", "designHR", "n.sim")

# Deleting rows/scenarios with equal censoring rate in control and treatment group
parameters <- parameters[-which(parameters$p.C.C==parameters$p.C.T),]

# Adding FU parameters, which are dependent on med.C
parameters$accrual.time <- 24 # 2 years = 24 months 
parameters$fu.time      <- 2*parameters$med.C

# Calculating sample size for each sub-scenario and adding n.sim to the 
# parameter matrix
parameters$n.obs.T <- NA
parameters$n.obs.C <- NA
for (i in 1:dim(parameters)[1]) {
  n <- func.sample.size(
    alpha        = alpha,
    beta         = parameters$beta[i],
    r            = parameters$r[i],
    designHR     = parameters$designHR[i],
    med.C        = parameters$med.C[i],
    fu.time      = parameters$fu.time[i],
    accrual.time = parameters$accrual.time[i],
    p.C          = parameters$p.C.C[i])
  parameters$n.obs.T[i] = n$n.obs.T
  parameters$n.obs.C[i] = n$n.obs.C
}

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path, "..parameters.rds"))


#### Data Generation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), .combine = rbind, 
                  .packages = c("survival")) %dopar% {
                    # Generating Data "n.sim" times for each sub-scenario
                    data.scen.i <- rep(list(NA), para$n.sim)
                    for (i in 1:para$n.sim) {
                      data.scen.i[[i]] <- 
                        DataEXP(
                          n.obs.T        = para$n.obs.T, 
                          n.obs.C        = para$n.obs.C,
                          median.control = para$med.C,
                          accrual.time   = para$accrual.time,
                          fu.time        = para$fu.time,
                          designHR       = para$designHR, 
                          HR.var         = para$HR.var, 
                          p.C            = para$p.C,
                          p.C.C          = para$p.C.C,
                          p.C.T          = para$p.C.T,
                          index.seed     = seed[i])
                    }
                    
                    # Saving data set
                    saveRDS(data.scen.i,
                            file = paste0(path,
                                          "..beta.",     para$beta,
                                          ".p.C.",       para$p.C,
                                          ".p.C.C.",     para$p.C.C,
                                          ".p.C.T.",     para$p.C.T,
                                          ".r.",         para$r,
                                          ".med.C.",     para$med.C,
                                          ".shape.C.T.", para$shape.C.T,
                                          ".HR.var.",    para$HR.var,
                                          ".designHR.",  para$designHR,
                                          ".accrual.",   para$accrual.time,
                                          ".fu.",        para$fu.time,
                                          ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

(time.1 <- Sys.time() - start.time.1)
