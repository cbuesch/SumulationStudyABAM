################################################################################
####    Functions needed for analysis of generated data of Simulation 1     ####
################################################################################

#----------------------- Analysis of generated trials --------------------------
DataAnalysis.Trial <- function(df){
  # sig = 0: Not significant
  # sig = 1: Significant
  # sig = 2: significant but no censoring in one of the treatment groups
  
  
  # Performing Cox Regression
  cox.sum <- summary(coxph(Surv(Data.TandC, Data.status != 0) ~ Data.Z, data = df))
  # Calculating censoring rate of the trial
  cens.rate <- sum(df$Data.status==0)/dim(df)[1]
  cens.rate.T <- sum(df$Data.status[df$Data.Z=="Treatment"]==0)/sum(df$Data.Z=="Treatment")
  cens.rate.C <- sum(df$Data.status[df$Data.Z=="Control"]==0)/sum(df$Data.Z=="Control")
  
  # Checking if log rank test is significant
  if(cox.sum$waldtest[[3]]>0.05 | cox.sum$conf.int[3]>1){
    # If log rank test is not significant, no additional benefit assessment
    # methods are applied and hence no analysis is performed
    Data.out <- data.frame(cens.rate, cens.rate.T, cens.rate.C,
                           HR.point = NA, HR.CI.low = NA, HR.CI.up = NA,
                           median.C = NA, median.T = NA, median.gain = NA,
                           surv.gain = NA, surv.C.std.err = NA, surv.T.std.err = NA,
                           IQWiG_RR = NA, ModIQWiG_HR = NA,
                           ESMO_WithoutSurvGain = NA, ESMO = NA, 
                           ASCO = NA, 
                           sig = 0)
  }else{
    # If log rank test is significant, additional benefit assessment
    # methods are applied
    
    ## Calculation of HR, upper and lower CI
    HR.point  <- cox.sum$conf.int[1]
    HR.CI.low <- cox.sum$conf.int[3]
    HR.CI.up  <- cox.sum$conf.int[4]
    
    ## Calculation of median survival times
    # If the survival curve does not fall below 50%, e.g. due to large
    # treatment effects, the median survival time cannot be observed.
    # In this case a conservative approach was implemented, using the last
    # observed censoring or event time point of the survival curve instead.
    survtest <- summary(survfit(Surv(Data.TandC, Data.status != 0) ~ Data.Z,
                                data = df))$table
    
    # Median of the treatment group (median.T)
    if(is.na(survtest[2,7])){
      median.T <- max(df$Data.TandC[df$Data.Z=="Treatment"])
    }else{
      median.T <- survtest[2,7]
    }
    
    # Median of the control group (median.C)
    if(is.na(survtest[1,7])){
      median.C <- max(df$Data.TandC[df$Data.Z=="Control"])
    }else{
      median.C <- survtest[1,7]
    }
    
    # Categorizing median.C --> decides which part of the ESMO is going to
    # be used
    if (median.C <=12){
      med.C = "<=12"
    }else if (median.C > 12 && median.C <= 24) {
      med.C = "12.24"
    } else if (median.C > 24) {
      med.C = ">24"
    }
    
    ## Calculation of median survival gain (median.gain)
    median.gain <- abs(median.T - median.C)
    
    ## Calculation of survival gain and standard error of bith treatment groups
    fit.C <- survfit(Surv(Data.TandC, Data.status != 0) ~ 1, data = df[which(df$Data.Z=="Control"),])
    fit.T <- survfit(Surv(Data.TandC, Data.status != 0) ~ 1, data = df[which(df$Data.Z=="Treatment"),])
    
    surv.gain <-
      switch (as.character(med.C),
              "<=12" = {
                surv.C = summary(fit.C, times=2*12, extend=TRUE)$surv
                surv.T = summary(fit.T, times=2*12, extend=TRUE)$surv
                if(is.null(surv.C) && is.null(surv.T)){
                  0
                }else if(is.null(surv.C)){
                  abs(0 - surv.T)
                }else if(is.null(surv.T)){
                  abs(surv.C - 0)
                }else{
                  abs(surv.C - surv.T)
                }
              },
              "12.24" = {
                surv.C = summary(fit.C, times=3*12, extend=TRUE)$surv
                surv.T = summary(fit.T, times=3*12, extend=TRUE)$surv
                if(is.null(surv.C) && is.null(surv.T)){
                  0
                }else if(is.null(surv.C)){
                  abs(0 - surv.T)
                }else if(is.null(surv.T)){
                  abs(surv.C - 0)
                }else{
                  abs(surv.C - surv.T)
                }
              },
              ">24" = {
                surv.C = summary(fit.C, times=5*12, extend=TRUE)$surv
                surv.T = summary(fit.T, times=5*12, extend=TRUE)$surv
                if(is.null(surv.C) && is.null(surv.T)){
                  0
                }else if(is.null(surv.C)){
                  abs(0 - surv.T)
                }else if(is.null(surv.T)){
                  abs(surv.C - 0)
                }else{
                  abs(surv.C - surv.T)
                }
              })
    
    surv.C.std.err <-
      switch (as.character(med.C),
              "<=12" = {
                summary(fit.C, times=2*12, extend=TRUE)$std.err
                # If survival rate is 0%, then std.error is NaN
              },
              "12.24" = {
                summary(fit.C, times=3*12, extend=TRUE)$std.err
                # If survival rate is 0%, then std.error is NaN
              },
              ">24" = {
                summary(fit.C, times=5*12, extend=TRUE)$std.err
                # If survival rate is 0%, then std.error is NaN
              })
    
    surv.T.std.err <-
      switch (as.character(med.C),
              "<=12" = {
                summary(fit.T, times=2*12, extend=TRUE)$std.err
                # If survival rate is 0%, then std.error is NaN
              },
              "12.24" = {
                summary(fit.T, times=3*12, extend=TRUE)$std.err
                # If survival rate is 0%, then std.error is NaN
              },
              ">24" = {
                summary(fit.T, times=5*12, extend=TRUE)$std.err
                # If survival rate is 0%, then std.error is NaN
              })
    
    ## Applying Additional Benefit Assessment Methods
    
    # IQWiG_RR
    IQWiG_RR <- (HR.CI.up<0.85)*6 + (HR.CI.up<0.95 & HR.CI.up>=0.85)*5 + (HR.CI.up<1 & HR.CI.up>=0.95)*4
    
    # Mod-IQWiG_HR
    ModIQWiG_HR <- (HR.CI.up<0.7908765)*6 + (HR.CI.up<0.9286669 & HR.CI.up>=0.7908765)*5 + (HR.CI.up<1 & HR.CI.up>=0.9286669)*4
    
    # ESMO
    ESMO_WithoutSurvGain <-
      switch (as.character(med.C),
              "<=12" = {
                (HR.CI.low>0.7 || median.gain<1.5)*1 +
                  (  (HR.CI.low<=0.65 && median.gain<2 && median.gain>=1.5) ||
                       (HR.CI.low<0.7 && HR.CI.low>0.65 && median.gain>=1.5)  )*2 +
                  (HR.CI.low<=0.65 && median.gain<3 && median.gain>=2)*3 +
                  (HR.CI.low<=0.65 && median.gain>=3)*4
              },
              "12.24" = {
                (HR.CI.low>0.75 || median.gain<1.5)*1 +
                  ( (HR.CI.low<=0.7 && median.gain<3 && median.gain>=1.5) ||
                      (HR.CI.low<0.75 && HR.CI.low>0.7 && median.gain>=1.5)  )*2 +
                  (HR.CI.low<=0.7 && median.gain<5 && median.gain>=3)*3 +
                  (HR.CI.low<=0.7 && median.gain>=5)*4
              },
              ">24" = {
                (HR.CI.low>0.75 || median.gain<4)*1 +
                  ( (HR.CI.low<=0.7 && median.gain<6 && median.gain>4) ||
                      (HR.CI.low<0.75 && HR.CI.low>0.7 && median.gain>=4)  )*2 +
                  (HR.CI.low<=0.7 && median.gain<9 && median.gain>=6)*3 +
                  (HR.CI.low<=0.7 && median.gain>=9)*4
              })
    
    ESMO <-
      switch (as.character(med.C),
              "<=12" = {
                if(surv.gain>=0.1){
                  4
                }else{
                  (HR.CI.low>0.7 || median.gain<1.5)*1 +
                    (  (HR.CI.low<=0.65 && median.gain<2 && median.gain>=1.5) ||
                         (HR.CI.low<0.7 && HR.CI.low>0.65 && median.gain>=1.5)  )*2 +
                    (HR.CI.low<=0.65 && median.gain<3 && median.gain>=2)*3 +
                    (HR.CI.low<=0.65 && median.gain>=3)*4
                }
              },
              "12.24" = {
                if(surv.gain>=0.1){
                  4
                }else{
                  (HR.CI.low>0.75 || median.gain<1.5)*1 +
                    ( (HR.CI.low<=0.7 && median.gain<3 && median.gain>=1.5) ||
                        (HR.CI.low<0.75 && HR.CI.low>0.7 && median.gain>=1.5)  )*2 +
                    (HR.CI.low<=0.7 && median.gain<5 && median.gain>=3)*3 +
                    (HR.CI.low<=0.7 && median.gain>=5)*4
                }
              },
              ">24" = {
                if(surv.gain>=0.1){
                  4
                }else{
                  (HR.CI.low>0.75 || median.gain<4)*1 +
                    ( (HR.CI.low<=0.7 && median.gain<6 && median.gain>4) ||
                        (HR.CI.low<0.75 && HR.CI.low>0.7 && median.gain>=4)  )*2 +
                    (HR.CI.low<=0.7 && median.gain<9 && median.gain>=6)*3 +
                    (HR.CI.low<=0.7 && median.gain>=9)*4
                }
              })
    
    # ASCO
    # HR Score (death)
    ASCO_HRScore <- (1 - HR.point)*100
    
    # TAIL OF THE CURVE
    # time point where survival curve is 2X the median OS of the comparator regimen (control)
    time_point <- median.C*2
    # proportion of patients alive at time_point
    help <- summary(survfit(Surv(Data.TandC, Data.status != 0) ~ Data.Z, data = df), times = time_point, extend = TRUE)
    p_C  <- help$surv[1]
    p_T  <- help$surv[2]
    if( (p_C>=0.2) && # assuming > 20% surviving with standard, otherwise no bonus points
        ( (p_T - p_C) / p_C >= 1.5) # 50% or greater improvement in proportion of patients alive with the test regimen
    ){
      ASCO <- ASCO_HRScore + 20
    }else{
      ASCO <- ASCO_HRScore
    }
    
    ## Return analysed data
    if (sum(df$Data.status[df$Data.Z=="Treatment"]==0)==0 ||
        sum(df$Data.status[df$Data.Z=="Control"]==0)==0) {
      # Significant but no censoring in one of the treatment groups
      Data.out <- data.frame(cens.rate, cens.rate.T, cens.rate.C,
                             HR.point, HR.CI.low, HR.CI.up,
                             median.C, median.T, median.gain,
                             surv.gain, surv.C.std.err, surv.T.std.err,
                             IQWiG_RR, ModIQWiG_HR,
                             ESMO_WithoutSurvGain, ESMO,
                             ASCO, 
                             sig = 2)
    } else{
      Data.out <- data.frame(cens.rate, cens.rate.T, cens.rate.C,
                             HR.point, HR.CI.low, HR.CI.up,
                             median.C, median.T, median.gain,
                             surv.gain, surv.C.std.err, surv.T.std.err,
                             IQWiG_RR, ModIQWiG_HR,
                             ESMO_WithoutSurvGain, ESMO,
                             ASCO,
                             sig = 1)
    }
  }
  return(Data.out)
}

#----------------------- Combining analysis results  ---------------------------
ScenarioCombining.NSim <- function(path_load_ana, para){
  
  ## Loading analysed trials of sub-scenario
  df.raw  <- readRDS(paste0(path_load_ana,
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
  
  ## Calculating power
  pow <- sum(df.raw$sig%in% c(1,2))/para$n.sim
  # Significant trials with censoring in both of the treatment groups
  pow_WithCensoringBothGroup <- sum(df.raw$sig==1)/para$n.sim
  # Significant trials without censoring in one of the treatment groups
  pow_WithoutCensoringBothGroup <- sum(df.raw$sig==2)/para$n.sim
  
  ## Selecting only significant trials for further analysis
  index.sig <- which(df.raw$sig %in% c(1,2))
  df.raw.help <- df.raw[index.sig,]
  
  ## Selecting only trials where survival median of control and treatment group can be calculated
  index.median <- which(!is.na(df.raw.help$median.gain))
  df.raw.sig <- df.raw.help[index.median,]
  # Survival median of Control group cannot be calculated
  perc.median.C.NA <- sum(is.na(df.raw.help$median.C))/dim(df.raw.help)[1]
  # Survival median of Treatment group cannot be calculated
  perc.median.T.NA <- sum(is.na(df.raw.help$median.T))/dim(df.raw.help)[1]
  
  ## Calculating the mean of continuous measures over all n.sim
  mean.CensRate     <- mean(df.raw.sig$cens.rate)
  mean.CensRate.T   <- mean(df.raw.sig$cens.rate.T)
  mean.CensRate.C   <- mean(df.raw.sig$cens.rate.C)
  mean.HRpoint      <- mean(df.raw.sig$HR.point)
  mean.HRCIlow      <- mean(df.raw.sig$HR.CI.low)
  mean.HRCIup       <- mean(df.raw.sig$HR.CI.up)
  mean.HRCIWidth    <- mean(df.raw.sig$HR.CI.up - df.raw.sig$HR.CI.low)
  mean.logHRCIlow   <- mean(log(df.raw.sig$HR.CI.low))
  mean.logHRCIup    <- mean(log(df.raw.sig$HR.CI.up))
  mean.logHRCIWidth <- mean(log(df.raw.sig$HR.CI.up) - log(df.raw.sig$HR.CI.low))
  mean.MedianC      <- mean(df.raw.sig$median.C)
  mean.MedianT      <- mean(df.raw.sig$median.T)
  mean.MedianGain   <- mean(df.raw.sig$median.gain)
  mean.SurvGain     <- mean(df.raw.sig$surv.gain)
  mean.ASCO         <- mean(df.raw.sig$ASCO)

  ## Survival rate std.err
  # Mean of Kaplan-Meier standard error at 2/3/5 years (with extensions)
  mean.Surv.C.std.err.extend                    <- mean(df.raw.sig$surv.C.std.err, na.rm = TRUE)
  mean.Surv.C.std.err.extend.IncludingSurvRate0 <- mean(ifelse(is.nan(df.raw.sig$surv.C.std.err), 0, df.raw.sig$surv.C.std.err))
  
  mean.Surv.T.std.err.extend                    <- mean(df.raw.sig$surv.T.std.err, na.rm = TRUE)
  mean.Surv.T.std.err.extend.IncludingSurvRate0 <- mean(ifelse(is.nan(df.raw.sig$surv.T.std.err), 0, df.raw.sig$surv.T.std.err))
  
  ## Maximal observed HR (point estimate) of significant trials
  max.HRpoint <- max(df.raw.sig$HR.point)
  
  ## IQWiG_RR
  # Percentage maximal IQWiG_RR grade
  df.raw.sig$IQWiG_RR.max <- df.raw.sig$IQWiG_RR==6
  perc.IQWiG_RR.max       <- sum(df.raw.sig$IQWiG_RR.max==1)/dim(df.raw.sig)[1]
  
  # Number of maximal and not maximal IQWiG_RR grade
  n.IQWiG_RR.max          <- sum(df.raw.sig$IQWiG_RR.max==1)
  n.IQWiG_RR.not.max      <- sum(df.raw.sig$IQWiG_RR.max!=1)
  
  # Calculation of maximal HR and corresponding CI for the maximum IQWiG_RR grade
  if(sum(df.raw.sig$IQWiG_RR.max) != 0){
    max.HRpoint.IQWiG_RR <- max(df.raw.sig[which(df.raw.sig$IQWiG_RR.max==1),]$HR.point)
    max.HRCIup.IQWiG_RR  <- max(df.raw.sig[which(df.raw.sig$IQWiG_RR.max==1),]$HR.CI.up)
    max.HRCIlow.IQWiG_RR <- max(df.raw.sig[which(df.raw.sig$IQWiG_RR.max==1),]$HR.CI.low)
  }else{
    max.HRpoint.IQWiG_RR <- NA
    max.HRCIup.IQWiG_RR  <- NA
    max.HRCIlow.IQWiG_RR <- NA
  }
  
  ## ModIQWiG_HR
  # Percentage maximal ModIQWiG_HR grade
  df.raw.sig$ModIQWiG_HR.max <- df.raw.sig$ModIQWiG_HR==6
  perc.ModIQWiG_HR.max       <- sum(df.raw.sig$ModIQWiG_HR.max==1)/dim(df.raw.sig)[1]
  
  # Number of maximal and not maximal ModIQWiG_HR grade
  n.ModIQWiG_HR.max   <- sum(df.raw.sig$ModIQWiG_HR.max==1)
  n.IQWIG.new.not.max <- sum(df.raw.sig$ModIQWiG_HR.max!=1)
  
  # Calculation of maximal HR and corresponding CI for the maximum ModIQWiG_HR grade
  if(sum(df.raw.sig$ModIQWiG_HR.max) != 0){
    max.HRpoint.IQWIG.new <- max(df.raw.sig[which(df.raw.sig$ModIQWiG_HR.max==1),]$HR.point)
    max.HRCIup.IQWIG.new  <- max(df.raw.sig[which(df.raw.sig$ModIQWiG_HR.max==1),]$HR.CI.up)
    max.HRCIlow.IQWIG.new <- max(df.raw.sig[which(df.raw.sig$ModIQWiG_HR.max==1),]$HR.CI.low)
  }else{
    max.HRpoint.IQWIG.new <- NA
    max.HRCIup.IQWIG.new  <- NA
    max.HRCIlow.IQWIG.new <- NA
  }
  
  ## ESMO without survival gain rule
  # Percentage of maximal grades
  perc.ESMO_WithoutSurvGain.max <- sum(df.raw.sig$ESMO_WithoutSurvGain==4)/dim(df.raw.sig)[1]
  
  # Number of maximal and not maximal maximal grade
  n.perc.ESMO_WithoutSurvGain.max     <-  sum(df.raw.sig$ESMO_WithoutSurvGain==4)
  n.perc.ESMO_WithoutSurvGain.not.max <-  sum(df.raw.sig$ESMO_WithoutSurvGain!=4)

  ## ESMO
  # Percentage of maximal grades due to surv.gain>=10%
  perc.ESMO.max.DueToSurvGain <- sum(df.raw.sig$surv.gain>=0.1)/sum(df.raw.sig$ESMO==4)
  
  # If maximal score achieved: why?
  ESMO.whyMax <- ifelse(
    df.raw.sig$ESMO_WithoutSurvGain==4 & df.raw.sig$surv.gain>=0.1, "Both",
    ifelse(
      df.raw.sig$ESMO_WithoutSurvGain!=4 & df.raw.sig$surv.gain>=0.1, "SurvGain",
      ifelse(
        df.raw.sig$ESMO_WithoutSurvGain==4 & df.raw.sig$surv.gain<0.1, "NotSurvGain",
        ifelse(
          df.raw.sig$ESMO_WithoutSurvGain!=4 & df.raw.sig$surv.gain<0.1, "None",
          NA
        )
      )
    )
  )
  perc.ESMO.whyMax.Both        <- sum(ESMO.whyMax=="Both")/sum(df.raw.sig$ESMO==4)
  perc.ESMO.whyMax.SurvGain    <- sum(ESMO.whyMax=="SurvGain")/sum(df.raw.sig$ESMO==4)
  perc.ESMO.whyMax.NotSurvGain <- sum(ESMO.whyMax=="NotSurvGain")/sum(df.raw.sig$ESMO==4)
  perc.ESMO.whyMax.None        <- sum(ESMO.whyMax=="None")/sum(df.raw.sig$ESMO==4)
  
  # Percentage maximal ESMO grade
  df.raw.sig$ESMO.max <- df.raw.sig$ESMO==4
  perc.ESMO.max       <- sum(df.raw.sig$ESMO==4)/dim(df.raw.sig)[1]
  
  # Number of maximal and not maximal ESMO grade
  n.perc.ESMO.max     <-  sum(df.raw.sig$ESMO==4)
  n.perc.ESMO.not.max <-  sum(df.raw.sig$ESMO!=4)
  
  # Calculation of maximal HR and corresponding CI for the maximum ESMO scores
  if(sum(df.raw.sig$ESMO.max) != 0){
    max.HRpoint.ESMO <- max(df.raw.sig[which(df.raw.sig$ESMO.max==1),]$HR.point)
    max.HRCIup.ESMO  <- max(df.raw.sig[which(df.raw.sig$ESMO.max==1),]$HR.CI.up)
    max.HRCIlow.ESMO <- max(df.raw.sig[which(df.raw.sig$ESMO.max==1),]$HR.CI.low)
  }else{
    max.HRpoint.ESMO <- NA
    max.HRCIup.ESMO  <- NA
    max.HRCIlow.ESMO <- NA
  }
  
  ## Percentages of maximal grades / significant trials compared to number of iterations
  perc.ModIQWiG_HR.max.n.sim <- sum(df.raw.sig$ModIQWiG_HR.max==1)/para$n.sim
  perc.IQWiG_RR.max.n.sim    <- sum(df.raw.sig$IQWiG_RR.max==1)/para$n.sim
  perc.ESMO.max.n.sim        <- sum(df.raw.sig$ESMO.max==1)/para$n.sim

  ## Pairwise correlation between additional benefit assessment methods
  ## (setting correlation to "n.a." if only the same grade was given form both
  ## methods which are compared)
  
  # ASCO / IQWiG_RR
  if(length(unique(df.raw.sig$IQWiG_RR)) == 1){
    Cor.ASCO.IQWiG_RR <- "n.a."
  }else{
    Cor.ASCO.IQWiG_RR <- cor(df.raw.sig$ASCO, df.raw.sig$IQWiG_RR, method = "spearman")
  }
  
  # ASCO / ModIQWiG_HR
  if(length(unique(df.raw.sig$ModIQWiG_HR)) == 1){
    Cor.ASCO.ModIQWiG_HR <- "n.a."
  }else{
    Cor.ASCO.ModIQWiG_HR <- cor(df.raw.sig$ASCO, df.raw.sig$ModIQWiG_HR, method = "spearman")
  }
  
  # ASCO / ESMO
  if(length(unique(df.raw.sig$ESMO)) == 1){
    Cor.ASCO.ESMO <- "n.a."
  }else{
    Cor.ASCO.ESMO  <- cor(df.raw.sig$ASCO, as.numeric(df.raw.sig$ESMO), method = "spearman")
  }
  
  # IQWiG_RR / ESMO
  if(length(unique(df.raw.sig$ESMO)) == 1 | length(unique(df.raw.sig$IQWiG_RR)) == 1){
    Cor.IQWiG_RR.ESMO <- "n.a."
  }else{
    Cor.IQWiG_RR.ESMO  <- cor(df.raw.sig$IQWiG_RR, as.numeric(df.raw.sig$ESMO), method = "spearman")
  }
  
  # ModIQWiG_HR / ESMO
  if(length(unique(df.raw.sig$ESMO)) == 1 | length(unique(df.raw.sig$ModIQWiG_HR)) == 1){
    Cor.ModIQWiG_HR.ESMO <- "n.a."
  }else{
    Cor.ModIQWiG_HR.ESMO  <- cor(df.raw.sig$ModIQWiG_HR, as.numeric(df.raw.sig$ESMO), method = "spearman")
  }
  
  # ModIQWiG_HR / IQWiG_RR
  if(length(unique(df.raw.sig$ModIQWiG_HR)) == 1 | length(unique(df.raw.sig$IQWiG_RR)) == 1 ){
    Cor.ModIQWiG_HR.IQWiG_RR <- "n.a."
  }else{
    Cor.ModIQWiG_HR.IQWiG_RR <- cor(df.raw.sig$ModIQWiG_HR, df.raw.sig$IQWiG_RR, method = "spearman")
  }
  
  ## Return combined analysed data
  results <- data.frame(pow, pow_WithCensoringBothGroup, pow_WithoutCensoringBothGroup,
                        
                        perc.median.C.NA, perc.median.T.NA,
                        
                        mean.CensRate, mean.CensRate.T, mean.CensRate.C,
                        
                        max.HRpoint, mean.HRpoint, mean.HRCIlow, mean.HRCIup,
                        mean.HRCIWidth, mean.logHRCIlow, mean.logHRCIup, mean.logHRCIWidth,
                        
                        mean.MedianC, mean.MedianT, mean.MedianGain, mean.SurvGain,
                        
                        mean.ASCO,
                        
                        mean.Surv.C.std.err.extend, mean.Surv.C.std.err.extend.IncludingSurvRate0,
                        mean.Surv.T.std.err.extend, mean.Surv.T.std.err.extend.IncludingSurvRate0,
                        
                        perc.IQWiG_RR.max, n.IQWiG_RR.max, n.IQWiG_RR.not.max,
                        max.HRpoint.IQWiG_RR, max.HRCIup.IQWiG_RR, max.HRCIlow.IQWiG_RR,
                        
                        perc.ModIQWiG_HR.max, n.ModIQWiG_HR.max, n.IQWIG.new.not.max,
                        max.HRpoint.IQWIG.new, max.HRCIup.IQWIG.new, max.HRCIlow.IQWIG.new,
                        
                        perc.ESMO_WithoutSurvGain.max, n.perc.ESMO_WithoutSurvGain.max, n.perc.ESMO_WithoutSurvGain.not.max,

                        perc.ESMO.max.DueToSurvGain,
                        perc.ESMO.whyMax.Both, perc.ESMO.whyMax.SurvGain,
                        perc.ESMO.whyMax.NotSurvGain, perc.ESMO.whyMax.None,
                        
                        perc.ESMO.max, n.perc.ESMO.max, n.perc.ESMO.not.max,
                        max.HRpoint.ESMO, max.HRCIup.ESMO, max.HRCIlow.ESMO,
                        
                        perc.ModIQWiG_HR.max.n.sim, perc.IQWiG_RR.max.n.sim,
                        perc.ESMO.max.n.sim,
                        
                        Cor.ASCO.IQWiG_RR, Cor.ASCO.ModIQWiG_HR, Cor.ASCO.ESMO,
                        Cor.IQWiG_RR.ESMO, Cor.ModIQWiG_HR.ESMO, Cor.ModIQWiG_HR.IQWiG_RR
  )
  return(results)
}

#------------------------ Analysis of each Scenario ----------------------------
#---------- (using "DataAnalysis.Trial" and "ScenarioCombining.NSim") ----------
func.DataAnalysis.Scenario <- function(path_load, path_save, Scenario){
  ## Loading parameters
  parameters <- readRDS(paste0(path_load, "..parameters.rds"))
  
  ## Sorting parameters (reduce computing time)
  parameters <- parameters[order(parameters$n.obs.C, decreasing = TRUE),]
  
  ## Analysis of each of the "n.sim" times simulated trials
  # Setting up parallel computing
  cl <- makeCluster(num.cl) 
  registerDoParallel(cl)
  
  # Starting parallel computing
  foreach(para = iter(parameters, by='row'),
          .packages = "survival", .combine = rbind, 
          .export = c("DataAnalysis.Trial")) %dopar% {
            # Load data of sub-scenarios with "n.sim" times replicates / trials
            df.loaded <- readRDS(paste0(path_load,
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
            
            # Analysis of each of the "n.sim" times simulated trials
            ana.data.scen.1           <- DataAnalysis.Trial(df = df.loaded[[1]])
            ana.data.scen.i           <- data.frame(matrix(NA, ncol = dim(ana.data.scen.1)[2], nrow = para$n.sim))
            ana.data.scen.i[1,]       <- ana.data.scen.1
            colnames(ana.data.scen.i) <- names(ana.data.scen.1)
            
            for (i in 2:para$n.sim) {
              ana.data.scen.i[i,] <- DataAnalysis.Trial(df = df.loaded[[i]])
            }
            
            # Saving of analysis results of each of the "n.sim" times simulated trials
            saveRDS(ana.data.scen.i,
                    file = paste0(path_save,
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
  
  ## Combining the analysed results of the "n.sim" trials of each sub-scenario
  # Resorting parameters by reloading parameters
  parameters  <- readRDS(paste0(path_load, "..parameters.rds"))
  
  # Setting up parallel computing
  cl <- makeCluster(num.cl) 
  registerDoParallel(cl)
  
  # Starting parallel computing
  results_raw <- foreach(para = iter(parameters, by='row'),
                         .combine = rbind,
                         .packages = c("stringr", "irr"), 
                         .export = c("ScenarioCombining.NSim")) %dopar% {
                           return(
                             ScenarioCombining.NSim(
                               path_load_ana = path_save,
                               para = para)
                           )
                         }
  
  # Stopping parallel computing 
  stopCluster(cl)
  
  # Saving the combined analysis of each sub-scenario
  results <- cbind(Scenario, parameters, results_raw)
  
  saveRDS(results,
          file = paste0(path_save, "..CombinedResults.rds")
  )
}
