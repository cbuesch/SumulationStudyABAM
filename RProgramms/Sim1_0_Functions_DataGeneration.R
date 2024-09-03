################################################################################
####     Costume functions needed for data generation of Simulation 1       ####
################################################################################


#-------------------------- Sample Size Calculation ----------------------------
func.sample.size <- function(alpha, beta, r, designHR, med.C, fu.time, 
                             accrual.time, p.C) {
  # Calculation of sample size for survival data (Log-Rank Test) via Schoenfeld 
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
  # Events
  d.raw = (1+r)^2/r * (qnorm(1-alpha/2) + qnorm(1-beta))^2/(log(designHR)^2)
  # Rounding events up so it can be splitted by r
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
DataEXP = function(n.obs.T, n.obs.C, median.control, accrual.time, fu.time, 
                   designHR, HR.var, p.C, p.C.C, p.C.T, index.seed) {
  ### Parameter descriptions
  # n.obs.T:        number of patients in the treatment group
  # n.obs.C:        number of patients in the control group
  # median.control: median survival time of the control group
  # designHR:       hazard ratio of treatment against control 
  # HR.var:         Variation of the trueHR and designHR. 0 --> designHR = trueHR
  
  # p.C:            censoring rate 
  #                 (if the treatment groups have different censoring rates, 
  #                 set to "NA")
  # p.C.C:          censoring rate of the control group (if both treatment 
  #                 groups have the same censoring rates, set to "NA")
  # p.C.T:          censoring rate of the treatment group (if both treatment 
  #                 groups have the same censoring rates, set to "NA")
  # accrual.time:   Accrual time (if no administrative censoring wants to be 
  #                 generated, set to "NA")
  # fu.time:        Follow-Up time (if no administrative censoring wants to be 
  #                 generated, set to "NA")
  # index.seed:     seed for the simulation
  
  
  ### Data generation
  set.seed(index.seed)
  
  ## Event times
  mytimes.C <- rexp(n.obs.C, rate = log(2)/median.control)
  mytimes.T <- rexp(n.obs.T, rate = (designHR*HR.var)*log(2)/median.control)
  mytimes   <- c(mytimes.C, mytimes.T)
  
  
  ## Combination of administrative and exponential censoring 
  if( (!is.na(p.C) && !is.na(fu.time) && !is.na(accrual.time)) |
      (is.na(p.C) && !is.na(p.C.C) && !is.na(p.C.T) && !is.na(fu.time) && !is.na(accrual.time)) ){
    
    # Both treatment groups have same censoring proportion
    if(!is.na(p.C)){
      # Administrative censoring
      cens.times.AC <- runif(n.obs.C+n.obs.T, min = 0, max = accrual.time) + fu.time
      
      # Data frame
      Data.help <- data.frame(id = 1:(n.obs.C + n.obs.T), T = mytimes, Z = factor(c(rep("Control",n.obs.C), rep("Treatment",n.obs.T))),
                              C = cens.times.AC)
      Data.help$status <- factor(as.numeric(Data.help$T <= Data.help$C))
      Data.help$TandC  <- pmin(Data.help$T,Data.help$C)
      
      # Censoring rate of the remaining events "mytimes.AC"
      mytimes.AC    <- Data.help$T[which(Data.help$status==1)]
      # Calculating of needed censoring rate to receive overall censoring rate of p.C 
      current.CensoringRate <- 1 # of Data.help$status==1 --> always 100%
      needed.CensoringRate  <- (p.C * dim(Data.help)[1] - sum(Data.help$status==0)*current.CensoringRate) / sum(Data.help$status==1)
      # if administrative censoring rate is already larger than wanted censoring rate of p.C, no specific censoring rate is added
      if (needed.CensoringRate > 0){
        cens.times.Exp <- sapply(mytimes.AC, function(mytimes.AC, needed.CensoringRate){rexp(1, rate = -log(1-needed.CensoringRate)/mytimes.AC)}, needed.CensoringRate)
      }else{
        # if administrative censoring rate is already larger than wanted censoring rate of p.C, 
        # no specific censoring rate is added
        # --> Setting censoring added censoring rates to Infinity --> event rates are always smaller than specific censoring rates
        cens.times.Exp <- rep(Inf, length(mytimes.AC))
      }
      
      # Final data Frame
      Data <- Data.help
      Data$C[which(Data$status==1)] <- cens.times.Exp
      Data$status <- factor(as.numeric(Data$T <= Data$C))
      Data$TandC  <- pmin(Data$T,Data$C)
      
      Data.out <- data.frame(Data$TandC, Data$status, Data$Z)
      
      
      # Treatment groups have different censoring proportions
    }else if(is.na(p.C) && !is.na(p.C.C) && !is.na(p.C.T)){
      # Administrative censoring
      cens.times.C.AC <- runif(n.obs.C, min = 0, max = accrual.time) + fu.time
      cens.times.T.AC <- runif(n.obs.T, min = 0, max = accrual.time) + fu.time
      
      # Data frame
      Data.help.C <- data.frame(id = 1:n.obs.C, T = mytimes.C, Z = rep("Control",n.obs.C),
                                C = cens.times.C.AC)
      Data.help.C$status <- factor(as.numeric(Data.help.C$T <= Data.help.C$C))
      Data.help.C$TandC  <- pmin(Data.help.C$T,Data.help.C$C)
      
      Data.help.T <- data.frame(id = 1:n.obs.T, T = mytimes.T, Z = rep("Treatment",n.obs.T),
                                C = cens.times.T.AC)
      Data.help.T$status <- factor(as.numeric(Data.help.T$T <= Data.help.T$C))
      Data.help.T$TandC  <- pmin(Data.help.T$T,Data.help.T$C)
      
      # Censoring rate of the remaining events "mytimes.C.AC" and "mytimes.T.AC"
      mytimes.C.AC    <- Data.help.C$T[which(Data.help.C$status==1)]
      mytimes.T.AC    <- Data.help.T$T[which(Data.help.T$status==1)]
      
      if(length(mytimes.C.AC)!=0) { # Checking if any specific censoring is still needed
        # calculating needed censoring rate to receive overall censoring rate of p.C 
        current.CensoringRate.C <- 1 # of Data.help.C$status==1 --> always 100%
        needed.CensoringRate.C  <- (p.C.C * dim(Data.help.C)[1] - sum(Data.help.C$status==0)*current.CensoringRate.C) / sum(Data.help.C$status==1)
        # if administrative censoring rate is already larger than wanted censoring rate of p.C, no specific censoring rate is added
        if (needed.CensoringRate.C > 0){
          cens.times.C.Exp <- sapply(mytimes.C.AC, function(mytimes.C.AC, needed.CensoringRate.C){rexp(1, rate = -log(1-needed.CensoringRate.C)/mytimes.C.AC)}, needed.CensoringRate.C)
        }else{
          # if administrative censoring rate is already larger than wanted censoring rate of p.C, 
          # no specific censoring rate is added
          # --> Setting censoring added censoring rates to Infinity --> event rates are always smaller than specific censoring rates
          cens.times.C.Exp <- rep(Inf, length(mytimes.C.AC))
        }
      }else{cens.times.C.Exp<-NA}
      
      if(length(mytimes.T.AC)!=0) { # Checking if any specific censoring is still needed
        # calculating needed censoring rate to receive overall censoring rate of p.C (e.g. 60%)
        current.CensoringRate.T <- 1 # of Data.help.T$status==1 --> always 100%
        needed.CensoringRate.T  <- (p.C.T * dim(Data.help.T)[1] - sum(Data.help.T$status==0)*current.CensoringRate.T) / sum(Data.help.T$status==1)
        # if administrative censoring rate is already larger than wanted censoring rate of p.C, no specific censoring rate is added
        if (needed.CensoringRate.T > 0){
          cens.times.T.Exp <- sapply(mytimes.T.AC, function(mytimes.T.AC, needed.CensoringRate.T){rexp(1, rate = -log(1-needed.CensoringRate.T)/mytimes.T.AC)}, needed.CensoringRate.T)
        }else{
          # if administrative censoring rate is already larger than wanted censoring rate of p.C, 
          # no specific censoring rate is added
          # --> Setting censoring added censoring rates to Infinity --> event rates are always smaller than specific censoring rates
          cens.times.T.Exp <- rep(Inf, length(mytimes.T.AC))
        }
      }else{cens.times.T.Exp = NA}
      
      # Final data frame
      Data.C <- Data.help.C
      Data.C$C[which(Data.C$status==1)] <- cens.times.C.Exp
      Data.C$status <- factor(as.numeric(Data.C$T <= Data.C$C))
      Data.C$TandC  <- pmin(Data.C$T,Data.C$C)
      
      Data.T <- Data.help.T
      Data.T$C[which(Data.T$status==1)] <- cens.times.T.Exp
      Data.T$status <- factor(as.numeric(Data.T$T <= Data.T$C))
      Data.T$TandC  <- pmin(Data.T$T,Data.T$C)
      
      Data     <- rbind(Data.C, Data.T)
      Data.out <- data.frame(Data$TandC, Data$status, Data$Z)
    }
    
    ## Only exponential censoring 
  }else if(is.na(fu.time) && is.na(accrual.time)){
    # Both treatment groups have same censoring proportion
    if(!is.na(p.C)){
      cens.times <- sapply(mytimes, function(mytimes, p.C){rexp(1, rate = -log(1-p.C)/mytimes)}, p.C)
      
      # Treatment groups have different censoring proportions
    }else if(is.na(p.C) && !is.na(p.C.C) && !is.na(p.C.T)){
      cens.times.C <- sapply(mytimes.C, function(mytimes.C, p.C.C){rexp(1, rate = -log(1-p.C.C)/mytimes.C)}, p.C.C)
      cens.times.T <- sapply(mytimes.T, function(mytimes.T, p.C.T){rexp(1, rate = -log(1-p.C.T)/mytimes.T)}, p.C.T)  
      cens.times   <- c(cens.times.C, cens.times.T)
    }
    
    # Final data frame
    Data <- data.frame(id = 1:(n.obs.C + n.obs.T), T = mytimes, Z = factor(c(rep("Control",n.obs.C), rep("Treatment",n.obs.T))),
                       C = cens.times)
    Data$status <- factor(as.numeric(Data$T <= Data$C))
    Data$TandC  <- pmin(Data$T,Data$C)
    
    Data.out <- data.frame(Data$TandC, Data$status, Data$Z)
  }
  
  return(Data.out)
}

#------------- Data Generation of WEIBULL distributed failure times ------------
#--------------- with a combination of administrative censoring ----------------
#-------------------------- and exponential censoring --------------------------
DataWEIB <- function(n.obs.T, n.obs.C, median.control, accrual.time, fu.time, 
                     designHR, HR.var, p.C, shape.C.T, index.seed) {
  ### Parameter descriptions
  # n.obs.T:        Number of patients in treatment group
  # n.obs.C:        Number of patients in control group
  # median.control: Median survival time of the control group
  # designHR:       Hazard ratio of treatment against control 
  # HR.var:         Variation of the trueHR and designHR. 0 --> designHR = trueHR
  # shape.C.T:      Shape parameter for control and treatment group
  # p.C:            censoring rate, censoring exponential distributed
  # accrual.time:   Accrual time
  # fu.time:        Follow-Up time
  # index.seed:     Seed for the simulation
  
  
  ### Data generation
  set.seed(index.seed)
  
  ## Event times
  mytimes.C = rweibull(n.obs.C, 
                       shape = shape.C.T, 
                       scale = median.control/((log(2))^(1/shape.C.T)) )
  mytimes.T = rweibull(n.obs.C, 
                       shape = shape.C.T, 
                       scale = (median.control^shape.C.T / (designHR*HR.var*log(2)))^(1/shape.C.T) )
  mytimes = c(mytimes.C, mytimes.T)
  
  ## Combination of administrative and exponential censoring 
  ## (both treatment groups have same censoring proportion)
  # Administrative censoring
  cens.times.AC <- runif(n.obs.C+n.obs.T, min = 0, max = accrual.time) + fu.time
  
  # Help Data Frame
  Data.help <- data.frame(id = 1:(n.obs.C + n.obs.T), T = mytimes, Z = factor(c(rep("Control",n.obs.C), rep("Treatment",n.obs.T))),
                          C = cens.times.AC)
  Data.help$status <- factor(as.numeric(Data.help$T <= Data.help$C))
  Data.help$TandC  <- pmin(Data.help$T,Data.help$C)
  
  # Censoring rate of the remaining events "mytimes.AC"
  mytimes.AC    <- Data.help$T[which(Data.help$status==1)]
  # Calculating needed censoring rate to receive overall censoring rate of p.C 
  current.CensoringRate <- 1 # of Data.help$status==1 --> always 100%
  needed.CensoringRate  <- (p.C * dim(Data.help)[1] - sum(Data.help$status==0)*current.CensoringRate) / sum(Data.help$status==1)
  # if administrative censoring rate is already larger than wanted censoring rate of p.C, no specific censoring rate is added
  if (needed.CensoringRate > 0){
    cens.times.SC <- sapply(mytimes.AC, function(mytimes.AC, needed.CensoringRate){rexp(1, rate = -log(1-needed.CensoringRate)/mytimes.AC)}, needed.CensoringRate)
  }else{
    # if administrative censoring rate is already larger than wanted censoring rate of p.C, 
    # no specific censoring rate is added
    # --> Setting censoring added censoring rates to Infinity --> event rates are always smaller than specific censoring rates
    cens.times.SC <- rep(Inf, length(mytimes.AC))
  }
  
  ## Final data Frame
  Data <- Data.help
  Data$C[which(Data$status==1)] <- cens.times.SC
  Data$status <- factor(as.numeric(Data$T <= Data$C))
  Data$TandC  <- pmin(Data$T,Data$C)
  
  Data.out <- data.frame(Data$TandC, Data$status, Data$Z)
  
  return(Data.out)
}


#------------ Data Generation of GOMPERTZ distributed failure times ------------
#--------------- with a combination of administrative censoring ----------------
#------------------------- and exponential censoring ---------------------------
DataGOMP = function(n.obs.T, n.obs.C, median.control, accrual.time, fu.time,
                    designHR, HR.var, p.C, shape.C.T, index.seed) {
  ### Parameter descriptions
  # n.obs.T:        Number of patients in treatment group
  # n.obs.C:        Number of patients in control group
  # median.control: Median survival time of the control group
  # designHR:       Hazard ratio of treatment against control 
  # HR.var:         variation of the trueHR and designHR. 0 --> designHR = trueHR
  # shape.C.T:      Shape parameter for control and treatment group
  # p.C:            Censoring rate, censoring exponential distributed
  # accrual.time:   Accrual time
  # fu.time:        Follow-Up time
  # index.seed:     Seed for the simulation
  
  
  ### Data generation
  set.seed(index.seed)
  
  ## Event times
  mytimes.C = rgompertz(n.obs.C, 
                        shape = shape.C.T, 
                        rate = shape.C.T*log(2) / (exp(median.control*shape.C.T)-1) )
  mytimes.T = rgompertz(n.obs.C, 
                        shape = shape.C.T, 
                        rate = (designHR*HR.var)*shape.C.T*log(2) / (exp(median.control*shape.C.T)-1) )
  mytimes = c(mytimes.C, mytimes.T)
  
  ## Combination of administrative and exponential censoring 
  ## (both treatment groups have same censoring proportion)
  # Administrative censoring
  cens.times.AC <- runif(n.obs.C+n.obs.T, min = 0, max = accrual.time) + fu.time
  
  # Help Data Frame
  Data.help <- data.frame(id = 1:(n.obs.C + n.obs.T), T = mytimes, Z = factor(c(rep("Control",n.obs.C), rep("Treatment",n.obs.T))),
                          C = cens.times.AC)
  Data.help$status <- factor(as.numeric(Data.help$T <= Data.help$C))
  Data.help$TandC  <- pmin(Data.help$T,Data.help$C)
  
  # Censoring rate of the remaining events "mytimes.AC"
  mytimes.AC    <- Data.help$T[which(Data.help$status==1)]
  # calculating needed censoring rate to receive overall censoring rate of p.C 
  current.CensoringRate <- 1 # of Data.help$status==1 --> always 100%
  needed.CensoringRate  <- (p.C * dim(Data.help)[1] - sum(Data.help$status==0)*current.CensoringRate) / sum(Data.help$status==1)
  # if administrative censoring rate is already larger than wanted censoring rate of p.C, no specific censoring rate is added
  if (needed.CensoringRate > 0){
    cens.times.SC <- sapply(mytimes.AC, function(mytimes.AC, needed.CensoringRate){rexp(1, rate = -log(1-needed.CensoringRate)/mytimes.AC)}, needed.CensoringRate)
  }else{
    # if administrative censoring rate is already larger than wanted censoring rate of p.C, 
    # no specific censoring rate is added
    # --> Setting censoring added censoring rates to Infinity --> event rates are always smaller than specific censoring rates
    cens.times.SC <- rep(Inf, length(mytimes.AC))
  }
  
  ## Final data frame
  Data <- Data.help
  Data$C[which(Data$status==1)] <- cens.times.SC
  Data$status <- factor(as.numeric(Data$T <= Data$C))
  Data$TandC  <- pmin(Data$T,Data$C)
  
  Data.out = data.frame(Data$TandC, Data$status, Data$Z)
  return(Data.out)
}

#---------- Data Generation of exponential distributed failure times -----------
#-------------------------- (non-proportional hazards) -------------------------
#--------------- with a combination of administrative censoring ----------------
#------------------------- and exponential censoring ---------------------------
DataEXPNonProp = function(n.obs.T, n.obs.C, median.control, accrual.time, 
                          fu.time, designHR, HR.var, effect.start.T, p.C, 
                          index.seed) {
  ### Parameter descriptions
  # n.obs.T:        Number of patients in treatment group
  # n.obs.C:        Number of patients in control group
  # median.control: Median survival time of the control group
  # designHR:       Hazard ratio of treatment against control 
  # HR.var:         variation of the trueHR and designHR. 0 --> designHR = trueHR
  # effect.start.T: Start of treatment effect of treatment group 
  #                 (late treatment effect)
  # p.C:            Censoring rate, censoring exponential distributed
  # accrual.time:   Accrual time
  # fu.time:        Follow-Up time
  # index.seed:     Seed for the simulation
  
  
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
  # Administrative censoring
  cens.times.AC <- runif(n.obs.C+n.obs.T, min = 0, max = accrual.time) + fu.time
  
  # Data Frame
  Data.help <- data.frame(id = 1:(n.obs.C + n.obs.T), T = mytimes, Z = factor(c(rep("Control",n.obs.C), rep("Treatment",n.obs.T))),
                          C = cens.times.AC)
  Data.help$status <- factor(as.numeric(Data.help$T <= Data.help$C))
  Data.help$TandC  <- pmin(Data.help$T,Data.help$C)
  
  # Specific censonring rate of the remaining events "mytimes.AC"
  mytimes.AC    <- Data.help$T[which(Data.help$status==1)]
  # calculating needed censoring rate to recieve overall censoring rate of p.C 
  current.CensoringRate <- 1 # of Data.help$status==1 --> always 100%
  needed.CensoringRate  <- (p.C * dim(Data.help)[1] - sum(Data.help$status==0)*current.CensoringRate) / sum(Data.help$status==1)
  # if administartive censoring rate is already larger than wanted censoring rate of p.C, no specific censoring rate is added
  if (needed.CensoringRate > 0){
    cens.times.SC <- sapply(mytimes.AC, function(mytimes.AC, needed.CensoringRate){rexp(1, rate = -log(1-needed.CensoringRate)/mytimes.AC)}, needed.CensoringRate)
  }else{
    # if administartive censoring rate is already larger than wanted censoring rate of p.C, 
    # no specific censoring rate is added
    # --> Setting censoring added censoring rates to Infinity --> event rates are always smaller than specific censoring rates
    cens.times.SC <- rep(Inf, length(mytimes.AC))
  }
  
  # Final data frame
  Data <- Data.help
  Data$C[which(Data$status==1)] <- cens.times.SC
  Data$status <- factor(as.numeric(Data$T <= Data$C))
  Data$TandC  <- pmin(Data$T,Data$C)
  
  Data.out <- data.frame(Data$TandC, Data$status, Data$Z)
  
  return(Data.out)
}

