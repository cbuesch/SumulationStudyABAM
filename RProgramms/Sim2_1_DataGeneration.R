################################################################################
######                 Data generation of Simulation 2                    ######
################################################################################

# Setting working directory (where functions and Simulation.seed are saved)
setwd("") 

#---------------------- Needed Packages and Functions --------------------------
library(foreach)
library(doParallel)

source("Sim2_0_Functions_DataGeneration.R")

#-------------------- Same parameter in all scenarios  -------------------------
# Number of used cores for Simulations
num.cl <- 35

# Number of iterations for each sub-scenario
n.sim <- 10000 

# Seed for the simulation:
#set.seed(1234567890)
#seed <- sample(x = 1:1000000000, size = n.sim, replace = FALSE)
#save(seed, file = "Sim2_Simulation.seed.RData")
load("Sim2_Simulation.seed.RData")


# Path for saving generated data:
path_data <- "" # e.g. ".../Simulation2/Data/" 


#------------------------- Standard Scenario 1 ---------------------------------
#### Parameters
path <- paste0(path_data, "Scen1")

med.C        <- c(6, 12, 18, 24, 30) 
designHR     <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
p.C          <- 0.6
shape.C.T    <- NA

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, p.C, r, med.C, shape.C.T, HR.var, designHR, n.sim)
# Adding colnames to parameter matrix
colnames(parameters) <- c("beta", "p.C", "r", "med.C", "shape.C.T", "HR.var", 
                          "designHR", "n.sim")

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

## Calculating lambda.cens
parameters$lambda.cens <- apply(X = parameters, MARGIN = 1, FUN = function(x){
  rootSolve::uniroot.all(f = function(aMax, lambda1, p1, lambda2, p2, duration, p.C, lambdaC){
    res <- 1/aMax * (
      p1*lambda1/((lambda1+lambdaC)^2) * ( exp((aMax-duration)*(lambda1+lambdaC)) - exp(-duration*(lambda1+lambdaC)) ) +
        p2*lambda2/((lambda2+lambdaC)^2) * ( exp((aMax-duration)*(lambda2+lambdaC)) - exp(-duration*(lambda2+lambdaC)) )
    ) + 
      p1 + p2 - p1*lambda1/(lambda1+lambdaC) - p2*lambda2/(lambda2+lambdaC) - 
      p.C
    return(res)
  }, 
  interval  = c(-1,1), 
  aMax      = x["accrual.time"], 
  lambda1   = log(2)/x["med.C"],
  p1        = 0.5,
  lambda2   = (x["designHR"]*x["HR.var"])*log(2)/x["med.C"], 
  p2        = 0.5,
  duration  = x["accrual.time"] + x["fu.time"],
  p.C = x["p.C"])
})

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
                          lambda.cens    = para$lambda.cens,
                          index.seed     = seed[i])
                    }
                    
                    # Saving data set
                    saveRDS(data.scen.i,
                            file = paste0(path,
                                          "..beta.",     para$beta,
                                          ".p.C.",       para$p.C,
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
p.C          <- 0.6
shape.C.T    <- NA 

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, p.C, r, med.C, shape.C.T, HR.var, designHR, n.sim)

# Adding colnames to parameter matrix
colnames(parameters) <- c("beta", "p.C", "r", "med.C", "shape.C.T", "HR.var",
                          "designHR", "n.sim")

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

# Calculating lambda.cens
parameters$lambda.cens <- apply(X = parameters, MARGIN = 1, FUN = function(x){
  rootSolve::uniroot.all(f = function(aMax, lambda1, p1, lambda2, p2, duration, p.C, lambdaC){
    res <- 1/aMax * (
      p1*lambda1/((lambda1+lambdaC)^2) * ( exp((aMax-duration)*(lambda1+lambdaC)) - exp(-duration*(lambda1+lambdaC)) ) +
        p2*lambda2/((lambda2+lambdaC)^2) * ( exp((aMax-duration)*(lambda2+lambdaC)) - exp(-duration*(lambda2+lambdaC)) )
    ) + 
      p1 + p2 - p1*lambda1/(lambda1+lambdaC) - p2*lambda2/(lambda2+lambdaC) - 
      p.C
    return(res)
  }, 
  interval  = c(-1,1), 
  aMax      = x["accrual.time"], 
  lambda1   = log(2)/x["med.C"],
  p1        = 0.5,
  lambda2   = (x["designHR"]*x["HR.var"])*log(2)/x["med.C"], 
  p2        = 0.5,
  duration  = x["accrual.time"] + x["fu.time"],
  p.C = x["p.C"])
})

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
                          lambda.cens    = para$lambda.cens,
                          index.seed     = seed[i])
                    }
                    
                    # Saving data set
                    saveRDS(data.scen.i,
                            file = paste0(path,
                                          "..beta.",     para$beta,
                                          ".p.C.",       para$p.C,
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
p.C          <- 0.6
shape.C.T    <- c(0.5, 1.5) 

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, p.C, r, med.C, shape.C.T, HR.var, designHR, n.sim)
# Adding colnames to parameter matrix
colnames(parameters) <- c("beta", "p.C", "r", "med.C", 
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

# Calculating lambda.cens
help_fiewerintegrate <- function(aMax, dur, p1, k1, l1, p2, k2, l2, p.C, lC){
  non.ad.cens.2 <- 1/aMax * 
    integrate(function(a, dur, p1, k1, l1, p2, k2, l2, lC){
      sapply(a, 
             function(a){
               integrate(function(t, a, dur, p1, k1, l1, p2, k2, l2, lC) {
                 return(
                   p1*k1*l1*(t*l1)^(k1-1) * exp(-(t*l1)^k1) * exp(-lC*t) + 
                     p2*k2*l2*(t*l2)^(k2-1) * exp(-(t*l2)^k2) * exp(-lC*t)
                 )
               }, lower = 0, upper = dur-a,
               a = a, dur = dur, 
               p1 = p1, k1 = k1, l1 = l1, 
               p2 = p2, k2 = k2, l2 = l2, lC = lC)$value
             }
      )
    }, lower = 0, upper = aMax,
    dur = dur, 
    p1 = p1, k1 = k1, l1 = l1, 
    p2 = p2, k2 = k2, l2 = l2, lC = lC)$value
  
  return(p1 + p2 - non.ad.cens.2 - p.C)
}

help_onlyIntegrate <- function(aMax, dur, p1, k1, l1, p2, k2, l2, p.C, lC){
  ad.cens <- 1/aMax * integrate(function(a, dur, p1, k1, l1, p2, k2, l2) {
    return(
      p1*(exp(-((dur-a)*l1)^k1)) + p2*(exp(-((dur-a)*l2)^k2))
    )
  }, lower = 0, upper = aMax,
  dur = dur, 
  p1 = p1, k1 = k1, l1 = l1, 
  p2 = p2, k2 = k2, l2 = l2)$value
  
  
  non.ad.cens <- 1/aMax * integrate(function(a, dur, p1, k1, l1, p2, k2, l2, lC){
    sapply(a, 
           function(a){
             integrate(function(t, a, dur, p1, k1, l1, p2, k2, l2, lC) {
               return(
                 (p1*(k1*l1*(t*l1)^(k1-1) * exp(-(t*l1)^k1)) + 
                    p2*(k2*l2*(t*l2)^(k2-1) * exp(-(t*l2)^k2))) * (1-exp(-lC*t))
               )
             }, lower = 0, upper = dur-a,
             a = a, dur = dur, 
             p1 = p1, k1 = k1, l1 = l1, 
             p2 = p2, k2 = k2, l2 = l2, lC = lC)$value
           }
    )
  }, lower = 0, upper = aMax,
  dur = dur, 
  p1 = p1, k1 = k1, l1 = l1, 
  p2 = p2, k2 = k2, l2 = l2, lC = lC)$value
  
  return(non.ad.cens + ad.cens - p.C)
}


parameters$lambda.cens <- apply(X = parameters, MARGIN = 1, FUN = function(x){
  uniroot(f = help_fiewerintegrate, 
          interval  = c(0,1), extendInt = "yes",
          aMax      = x["accrual.time"], 
          dur       = x["accrual.time"] + x["fu.time"],
          p1        = 0.5,
          k1        = x["shape.C.T"],
          l1        = ((log(2))^(1/x["shape.C.T"])) / x["med.C"],
          p2        = 0.5,
          k2        = x["shape.C.T"],
          l2        = ((x["designHR"]*x["HR.var"]*log(2)) / (x["med.C"]^x["shape.C.T"]))^(1/x["shape.C.T"]), 
          p.C       = x["p.C"])$root
})

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
                          lambda.cens    = para$lambda.cens,
                          shape.C.T      = para$shape.C.T,
                          index.seed     = seed[i])
                    }
                    
                    # Saving data set
                    saveRDS(data.scen.i,
                            file = paste0(path,
                                          "..beta.",     para$beta,
                                          ".p.C.",       para$p.C,
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
p.C          <- 0.6
shape.C.T    <- c(-0.2, 0.2) 

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, p.C, r, med.C, shape.C.T, HR.var, designHR, n.sim)

# Adding colnames to parameter matrix
colnames(parameters) <- c("beta", "p.C", "r", "med.C", "shape.C.T", "HR.var",
                          "designHR", "n.sim")

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

# Calculating lambda.cens
help_fiewerIntegrate <- function(aMax, dur, p1, a1, b1, p2, a2, b2, p.C, lC){
  non.ad.cens2 <- 1/aMax * integrate(function(a, dur, p1, a1, b1, p2, a2, b2, lC){
    sapply(a, 
           function(a){
             integrate(function(t, a, dur, p1, a1, b1, p2, a2, b2, lC) {
               return(
                 (
                   p1*b1*exp(a1*t)*exp( -b1/a1*(exp(a1*t)-1) )  +
                     p2*b2*exp(a2*t)*exp( -b2/a2*(exp(a2*t)-1) )
                 ) * exp(-lC*t)
               )
             }, lower = 0, upper = dur-a,
             a = a, dur = dur, 
             p1 = p1, a1 = a1, b1 = b1, 
             p2 = p2, a2 = a2, b2 = b2, lC = lC)$value
           }
    )
  }, lower = 0, upper = aMax,
  dur = dur, 
  p1 = p1, a1 = a1, b1 = b1, 
  p2 = p2, a2 = a2, b2 = b2, lC = lC)$value
  
  return(p1+p2-non.ad.cens2 - p.C)
}

help_onlyIntegrate <- function(aMax, dur, p1, a1, b1, p2, a2, b2, p.C, lC){
  ad.cens <- 1/aMax * integrate(function(a, dur, p1, a1, b1, p2, a2, b2) {
    return(
      p1*( exp(-b1/a1* (exp(a1*(dur-a))-1) ) ) + 
        p2*( exp(-b2/a2* (exp(a2*(dur-a))-1) ) )
    )
  }, lower = 0, upper = aMax,
  dur = dur, 
  p1 = p1, a1 = a1, b1 = b1, 
  p2 = p2, a2 = a2, b2 = b2)$value
  
  
  non.ad.cens <- 1/aMax * integrate(function(a, dur, p1, a1, b1, p2, a2, b2, lC){
    sapply(a, 
           function(a){
             integrate(function(t, a, dur, p1, a1, b1, p2, a2, b2, lC) {
               return(
                 (
                   p1*( b1*exp(a1*t)*exp( -b1/a1*(exp(a1*t)-1) ) ) +
                     p2*( b2*exp(a2*t)*exp( -b2/a2*(exp(a2*t)-1) ) )
                 ) * (1-exp(-lC*t))
               )
             }, lower = 0, upper = dur-a,
             a = a, dur = dur, 
             p1 = p1, a1 = a1, b1 = b1, 
             p2 = p2, a2 = a2, b2 = b2, lC = lC)$value
           }
    )
  }, lower = 0, upper = aMax,
  dur = dur, 
  p1 = p1, a1 = a1, b1 = b1, 
  p2 = p2, a2 = a2, b2 = b2, lC = lC)$value
  
  return(non.ad.cens + ad.cens - p.C)
}


parameters$lambda.cens <- apply(X = parameters, MARGIN = 1, FUN = function(x){
  uniroot(f = help_fiewerIntegrate, 
          interval  = c(0,1), extendInt = "yes",
          aMax      = x["accrual.time"], 
          dur       = x["accrual.time"] + x["fu.time"],
          p1        = 0.5,
          a1        = x["shape.C.T"],
          b1        = x["shape.C.T"]*log(2) / (exp(x["med.C"]*x["shape.C.T"])-1),
          p2        = 0.5,
          a2        = x["shape.C.T"],
          b2        = (x["designHR"]*x["HR.var"])*x["shape.C.T"]*log(2) / (exp(x["med.C"]*x["shape.C.T"])-1),
          p.C       = x["p.C"])$root
})

# Excluding scenarios where lambda.cens < 0 (administrative censoring already 
# too strong so that a larger censoring rate than wanted is present):
parameters <- parameters[which(parameters$lambda.cens>0),]

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
                          lambda.cens    = para$lambda.cens,
                          shape.C.T      = para$shape.C.T,
                          index.seed     = seed[i])
                    }
                    
                    # Saving data set
                    saveRDS(data.scen.i,
                            file = paste0(path,
                                          "..beta.",     para$beta,
                                          ".p.C.",       para$p.C,
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
p.C          <- 0.6
shape.C.T    <- NA
# treat.effect.start is set to 1/3*med.C due to ESMOs design of using "gain" (med.T-med.C)

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, p.C, r, med.C, shape.C.T, HR.var, designHR, n.sim)

# Adding colnames to parameter matrix
colnames(parameters) <- c("beta", "p.C", "r", "med.C", "shape.C.T", "HR.var", 
                          "designHR", "n.sim")

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

#  Calculating lambda.cens.C and lambda.cens.T
help_Integrate <- function(aMax, dur, startT, p1, l1, p2, l2, p.C, lC){
  ad.cens <- 1/aMax * integrate(function(a, dur, startT, p1, l1, p2, l2) {
    res <- ifelse((dur-a) > startT,
                  p1*exp(-l1*(dur-a)) + p2*( exp(-l1*startT) * exp(-l2*((dur-a)-startT)) ),
                  p1*exp(-l1*(dur-a)) + p2*exp(-l2*(dur-a))
    )
    return(res)
  }, lower = 0, upper = aMax,
  dur = dur, startT = startT,
  p1 = p1, l1 = l1,
  p2 = p2, l2 = l2)$value
  
  
  non.ad.cens <- 1/aMax * integrate(function(a, dur, startT, p1, l1, p2, l2, lC){
    res <-
      sapply(a,
             function(a){
               integrate(function(t, a, dur, startT, p1, l1, p2, l2, lC) {
                 ifelse(t>startT,
                        ( p1*l1*exp(-l1*t) + p2*l2*exp(-l1*startT)*exp(-l2*(t-startT)) ) * (1-exp(-lC*t)),
                        ( p1*l1*exp(-l1*t) + p2*l2*exp(-l2*t) ) * (1-exp(-lC*t))
                 )
               }, lower = 0, upper = dur-a,
               a = a, dur = dur, startT = startT,
               p1 = p1, l1 = l1,
               p2 = p2, l2 = l2, lC = lC)$value
             }
      )
    return(res)
  }, lower = 0, upper = aMax,
  dur = dur, startT = startT,
  p1 = p1, l1 = l1,
  p2 = p2, l2 = l2, lC = lC)$value
  
  return(non.ad.cens + ad.cens - p.C)
}

help_NoIntegrate <- function(aMax, dur, startT, p1, l1, p2, l2, p.C, lC){
  cens <- 1/aMax * (
    p1*l1/((l1+lC)^2) * ( exp(-(dur-aMax)*(l1+lC)) - exp(-dur*(l1+lC)) )
    + p2*l2/((l2+lC)^2) * exp(-l1*startT) * ( exp(-l2*(dur-aMax-startT) - lC*(dur-aMax)) - exp(-l2*(dur-startT) - lC*dur) )
  ) + p1 - p2*( exp(-l2*startT) - exp(-l1*startT) - 1 ) - p1*l1/(l1+lC) - p2*l2/(l2+lC) * ( exp(-startT*(l1 + lC)) - exp(-startT*(l2+lC)) + 1 )
  
  return(cens - p.C)
}

parameters$lambda.cens <- apply(X = parameters, MARGIN = 1, FUN = function(x){
  uniroot(f = help_NoIntegrate, 
          interval  = c(0,1), extendInt = "yes",
          aMax      = x["accrual.time"], 
          dur       = x["accrual.time"] + x["fu.time"],
          startT    = 1/3 * x["med.C"],
          p1        = 0.5,
          l1        = log(2) / x["med.C"],
          p2        = 0.5,
          l2        = (x["designHR"]*x["HR.var"])*log(2) / x["med.C"],
          p.C       = x["p.C"])$root
})

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
        lambda.cens    = para$lambda.cens,
        index.seed     = seed[i])
  }
  
  # Saving data set
  saveRDS(data.scen.i,
          file = paste0(path,
                        "..beta.",     para$beta,
                        ".p.C.",       para$p.C,
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