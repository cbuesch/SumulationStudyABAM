################################################################################
####     Costume functions needed for data generation of Simulation 2       ####
################################################################################



#-------------------------- Sample Size Calculation ----------------------------
func.sample.size <- function(alpha, beta, r, designHR, med.C, fu.time, 
                             accrual.time, p.C) {
  # Calculation of sample size for survival data (Log-Rank Test) via Schonfeld 
  # approach
  
  ### Parameter descriptions
  # alpha:        Type-I-error
  # beta:         Type-II-error
  # r:            Sample size proportion of the two intervention groups: n.T/n.C
  # designHR:     Hazard ratio of treatment against control 
  # med.C:        Median survival time of the control group
  # fu.time:      Follow-up time 
  # accrual.time: Accrual time 
  # p.C:          Censoring rate, censoring exponential distributed
  
  ### Sample size calculation
  ## Lambdas for exponential distribution:   
  lambda.C = log(2)/med.C
  lambda.T = designHR*log(2)/med.C   
  
  
  ## Shoenefeld:
  # events
  d.raw = (1+r)^2/r * (qnorm(1-alpha/2) + qnorm(1-beta))^2/(log(designHR)^2)
  # rounding events up so it can be splitted by r
  d = ceiling(d.raw)-1
  response = 1
  while (response!=0) { 
    d  = d + 1
    response = d %% (r+1)
  }
  
  # Sample size
  # Combination of administrative and exponential censoring
  if(!is.na(p.C) && !is.na(fu.time) && !is.na(accrual.time)){
    P.D = 1 - 1/(6*(1+r)) * (
      exp(-lambda.C*fu.time) + r*exp(-lambda.T*fu.time) +
        4*( exp(-lambda.C*(accrual.time/2 + fu.time)) + r*exp(-lambda.T*(accrual.time/2 + fu.time)) ) + 
        exp(-lambda.C*(accrual.time+fu.time)) + r*exp(-lambda.T*(accrual.time+fu.time))
    )
    # if probability of event using administrative censoring is smaller than (1-p.C), 
    # than use only P.D for calculation of n.raw, otherwise use only p.C 
    # In other words: If censoring rate of administrative censoring is larger than from specific censoring, 
    #                 than use only p.C for calculation of n.raw, otherwise use only P.D 
    if(P.D < (1-p.C)){
      n.raw = d/P.D
    }else{
      n.raw = d/(1-p.C)
    }
    
    # Exponential censoring (without administrative censoring)
  }else if(!is.na(p.C)){
    n.raw = d/(1-p.C)
    
    # Administrative data generation
  }else if(!is.na(fu.time) && !is.na(accrual.time)){
    # Probability of an event
    P.D = 1 - 1/(6*(1+r)) * (
      exp(-lambda.C*fu.time) + r*exp(-lambda.T*fu.time) +
        4*( exp(-lambda.C*(accrual.time/2 + fu.time)) + r*exp(-lambda.T*(accrual.time/2 + fu.time)) ) + 
        exp(-lambda.C*(accrual.time+fu.time)) + r*exp(-lambda.T*(accrual.time+fu.time))
    )
    n.raw = d/P.D
  }
  
  # Rounding sample size up so it can be splitted by r
  n.schoen = ceiling(n.raw)-1
  response = 1
  while (response!=0) { 
    n.schoen  = n.schoen + 1
    response = n.schoen %% (r+1)
  }
  
  return(list(n.obs.T = r*n.schoen/(r+1), n.obs.C = n.schoen/(r+1)))
}

#----------- Data Generation of exponential distributed failure times ----------
#--------------- with a combination of administrative censoring ----------------
#-------------------- and independent exponential censoring --------------------
DataEXP = function(n.obs.T, n.obs.C, median.control, accrual.time, fu.time, 
                   designHR, HR.var, lambda.cens, index.seed) {
  # Generating of survival Data 
  # n.obs.T:        Number of patients in treatment group
  # n.obs.C:        Number of patients in control group
  # median.control: Median survival time of the control group
  # designHR:       Hazard ratio of treatment against control 
  # HR.var:         Variation of the trueHR and designHR. 0 --> designHR = trueHR
  # lambda.cens:    Lambda of exponential distributed censoring
  # accrual.time:   Accrual time
  # fu.time:        Follow-Up time 
  # index.seed:     seed for the simulation
  
  ### Data generation
  set.seed(index.seed)
  
  ## Event times
  mytimes.C <- rexp(n.obs.C, rate = log(2)/median.control)
  mytimes.T <- rexp(n.obs.T, rate = (designHR*HR.var)*log(2)/median.control)
  mytimes   <- c(mytimes.C, mytimes.T)
  
  ## Combination of administrative and exponential censoring 
  cens.time <- c(rexp(n.obs.C, rate = lambda.cens),
                 rexp(n.obs.T, rate = lambda.cens))
  acc.time  <- runif(n.obs.C+n.obs.T, min = 0, max = accrual.time)
  
  ## Final data frame
  Data <- data.frame(id = 1:(n.obs.C + n.obs.T), T = mytimes, 
                     Z = factor(c(rep("Control",n.obs.C), rep("Treatment",n.obs.T))),
                     cens.time = cens.time, acc.time = acc.time)
  Data$TandC  <- pmin(Data$T, Data$cens.time, fu.time+accrual.time - Data$acc.time)
  Data$status <- factor(ifelse(Data$T<=pmin(Data$cens.time, fu.time+accrual.time - Data$acc.time), 1, 0))
  
  Data.out <- data.frame(Data$TandC, Data$status, Data$Z)
  
  return(Data.out)
}

#------------- Data Generation of WEIBULL distributed failure times ------------
#--------------- with a combination of administrative censoring ----------------
#-------------------- and independent exponential censoring --------------------
DataWEIB = function(n.obs.T, n.obs.C, median.control, accrual.time, fu.time, 
                    designHR, HR.var, lambda.cens, shape.C.T, index.seed) {
  ### Parameter descriptions
  # n.obs.T:        Number of patients in treatment group
  # n.obs.C:        Number of patients in control group
  # median.control: Median survival time of the control group
  # designHR:       Hazard ratio of treatment against control 
  # HR.var:         Variation of the trueHR and designHR. 0 --> designHR = trueHR
  # shape.C.T:      Shape parameter for control and treatment group
  # lambda.cens:    Lambda of exponential distributed censoring
  # accrual.time:   Accrual time
  # fu.time:        Follow-Up time
  # index.seed:     seed for the simulation
  
  
  ### Data generation
  set.seed(index.seed)
  
  ## Event times
  mytimes.C = rweibull(n.obs.C, shape = shape.C.T, scale = median.control/((log(2))^(1/shape.C.T)) )
  mytimes.T = rweibull(n.obs.C, shape = shape.C.T, scale = (median.control^shape.C.T / (designHR*HR.var*log(2)))^(1/shape.C.T) )
  mytimes = c(mytimes.C, mytimes.T)
  
  ## Combination of administrative and exponential censoring 
  cens.time <- c(rexp(n.obs.C, rate = lambda.cens),
                 rexp(n.obs.T, rate = lambda.cens))
  acc.time  <- runif(n.obs.C+n.obs.T, min = 0, max = accrual.time)
  
  # Final data frame
  Data <- data.frame(id = 1:(n.obs.C + n.obs.T), T = mytimes, Z = factor(c(rep("Control",n.obs.C), rep("Treatment",n.obs.T))),
                     cens.time = cens.time, acc.time = acc.time)
  Data$TandC  <- pmin(Data$T, Data$cens.time, fu.time+accrual.time - Data$acc.time)
  Data$status <- factor(ifelse(Data$T<=pmin(Data$cens.time, fu.time+accrual.time - Data$acc.time), 1, 0))
  
  Data.out <- data.frame(Data$TandC, Data$status, Data$Z)
  
  return(Data.out)
}

#------------ Data Generation of GOMPERTZ distributed failure times ------------
#--------------- with a combination of administrative censoring ----------------
#-------------------- and independent exponential censoring --------------------
DataGOMP <- function(n.obs.T, n.obs.C, median.control, accrual.time, fu.time, 
                     designHR, HR.var, lambda.cens, shape.C.T, index.seed) {
  ### Parameter descriptions
  # n.obs.T:        Number of patients in treatment group
  # n.obs.C:        Number of patients in control group
  # median.control: Median survival time of the control group
  # designHR:             hazard ratio of treatment against control 
  # HR.var:         Variation of the trueHR and designHR. 0 --> designHR = trueHR
  # shape.C.T:      Shape parameter for control and treatment group
  # lambda.cens:    Lambda of exponential distributed censoring
  # accrual.time:   Accrual time
  # fu.time:        Follow-Up time
  # index.seed:     seed for the simulation
  
  ### Data generation
  set.seed(index.seed)
  
  ## Event times
  mytimes.C = rgompertz(n.obs.C, shape = shape.C.T, rate = shape.C.T*log(2) / (exp(median.control*shape.C.T)-1) )
  mytimes.T = rgompertz(n.obs.C, shape = shape.C.T, rate = (designHR*HR.var)*shape.C.T*log(2) / (exp(median.control*shape.C.T)-1) )
  mytimes = c(mytimes.C, mytimes.T)
  
  ## Combination of administrative and exponential censoring 
  cens.time <- c(rexp(n.obs.C, rate = lambda.cens),
                 rexp(n.obs.T, rate = lambda.cens))
  acc.time  <- runif(n.obs.C+n.obs.T, min = 0, max = accrual.time)
  
  # Data Frame
  Data <- data.frame(id = 1:(n.obs.C + n.obs.T), T = mytimes, Z = factor(c(rep("Control",n.obs.C), rep("Treatment",n.obs.T))),
                     cens.time = cens.time, acc.time = acc.time)
  Data$TandC  <- pmin(Data$T, Data$cens.time, fu.time+accrual.time - Data$acc.time)
  Data$status <- factor(ifelse(Data$T<=pmin(Data$cens.time, fu.time+accrual.time - Data$acc.time), 1, 0))
  
  Data.out <- data.frame(Data$TandC, Data$status, Data$Z)
  
  return(Data.out)
}

#---------- Data Generation of exponential distributed failure times -----------
#-------------------------- (non-proportional hazards) -------------------------
#--------------- with a combination of administrative censoring ----------------
#-------------------- and independent exponential censoring --------------------
DataEXPNonProp <- function(n.obs.T, n.obs.C, median.control, accrual.time, 
                           fu.time, designHR, HR.var, effect.start.T, 
                           lambda.cens, index.seed) {
  ### Parameter descriptions
  # n.obs.T:        Number of patients in treatment group
  # n.obs.C:        Number of patients in control group
  # median.control: Median survival time of the control group
  # designHR:       Hazard ratio of treatment against control 
  # HR.var:         variation of the trueHR and designHR. 0 --> designHR = trueHR
  # effect.start.T: Start of treatment effect of treatment group 
  #                 (late treatment effect)
  # lambda.cens:    Lambda of exponential distributed censoring
  # accrual.time:   Accrual time
  # fu.time:        Follow-Up time 
  # index.seed:     seed for the simulation
  
  
  ### Data generation
  set.seed(index.seed)
  
  ## Event times
  # Control group
  mytimes.C <- rexp(n.obs.C, rate = log(2)/median.control)
  # Treatment group using inversion method
  # Inverse of probability distribution F
  Inv.Probabilityfunction.T <- function(y, lambda.C, lambda.T, treat.start){
    return(
      ifelse(y <= 1-exp(-lambda.C*treat.start),
             -log(1-y)/lambda.C,
             -(log(1-y) + lambda.C*treat.start)/lambda.T + treat.start)
    )
  }
  U <- runif(n.obs.T, min = 0, max = 1)
  mytimes.T <- Inv.Probabilityfunction.T(y = U, 
                                         lambda.C = log(2)/median.control, 
                                         lambda.T = (designHR*HR.var)*log(2)/median.control, 
                                         treat.start = effect.start.T)
  # Combining event times of both groups
  mytimes   <- c(mytimes.C, mytimes.T)
  
  
  
  ## Combination of administrative and exponential censoring 
  ## (both treatment groups have same censoring proportion)
  cens.time <- c(rexp(n.obs.C, rate = lambda.cens),
                 rexp(n.obs.T, rate = lambda.cens))
  acc.time  <- runif(n.obs.C+n.obs.T, min = 0, max = accrual.time)
  
  # Final data frame
  Data <- data.frame(id = 1:(n.obs.C + n.obs.T), T = mytimes, Z = factor(c(rep("Control",n.obs.C), rep("Treatment",n.obs.T))),
                     cens.time = cens.time, acc.time = acc.time)
  Data$TandC  <- pmin(Data$T, Data$cens.time, fu.time+accrual.time - Data$acc.time)
  Data$status <- factor(ifelse(Data$T<=pmin(Data$cens.time, fu.time+accrual.time - Data$acc.time), 1, 0))
  
  Data.out <- data.frame(Data$TandC, Data$status, Data$Z)
  
  return(Data.out)
}

