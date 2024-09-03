################################################################################
#########       Calculating ASCO cutoff values of Simulation 1         #########
################################################################################

#---------------------------- Needed packages ----------------------------------
library(tidyverse)
library(stringr)
library(cutpointr)
library(vcd)
library(data.table)

## Number of used cores for Analysis
num.cl <- 10

#-------------------- Defining loading/saving paths ----------------------------
path_start <- "" # e.g. ".../Simulation1/" 

# Saving path
path_saving <- paste0(path_start, "Ana/")

# Data loading paths
path_load <- tibble(
  Scenario = factor(
    c("Scen1", "Scen2", "Scen3WEIB",  "Scen3GOMP", "Scen4", "Scen5",
      "Scen6", "Scen7")
  ),
  path_load_data = NA,
  path_load_ana  = NA)

# Standard Scenario 1
path_load$path_load_data[which(path_load$Scenario=="Scen1")] <- paste0(path_start, "Data/Scen1")
path_load$path_load_ana[which(path_load$Scenario=="Scen1")]  <- paste0(path_start, "Ana/Scen1")
# Scenario 2
path_load$path_load_data[which(path_load$Scenario=="Scen2")] <- paste0(path_start, "Data/Scen2")
path_load$path_load_ana[which(path_load$Scenario=="Scen2")]  <- paste0(path_start, "Ana/Scen2")
# Scenario 3
path_load$path_load_data[which(path_load$Scenario=="Scen3WEIB")] <- paste0(path_start, "Data/Scen3WEIB")
path_load$path_load_ana[which(path_load$Scenario=="Scen3WEIB")]  <- paste0(path_start, "Ana/Scen3WEIB")
path_load$path_load_data[which(path_load$Scenario=="Scen3GOMP")] <- paste0(path_start, "Data/Scen3GOMP")
path_load$path_load_ana[which(path_load$Scenario=="Scen3GOMP")]  <- paste0(path_start, "Ana/Scen3GOMP")
# Scenario 4
path_load$path_load_data[which(path_load$Scenario=="Scen4")] <- paste0(path_start, "Data/Scen4")
path_load$path_load_ana[which(path_load$Scenario=="Scen4")]  <- paste0(path_start, "Ana/Scen4")
# Scenario 5
path_load$path_load_data[which(path_load$Scenario=="Scen5")] <- paste0(path_start, "Data/Scen5")
path_load$path_load_ana[which(path_load$Scenario=="Scen5")]  <- paste0(path_start, "Ana/Scen5")
# Scenario 6
path_load$path_load_data[which(path_load$Scenario=="Scen6")] <- paste0(path_start, "Data/Scen6")
path_load$path_load_ana[which(path_load$Scenario=="Scen6")]  <- paste0(path_start, "Ana/Scen6")
# Scenario 7
path_load$path_load_data[which(path_load$Scenario=="Scen7")] <- paste0(path_start, "Data/Scen7")
path_load$path_load_ana[which(path_load$Scenario=="Scen7")]  <- paste0(path_start, "Ana/Scen7")

#----------- Combining results of Scenarios and saving -------------------------------------
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)


# starting parallel computing
df.final = foreach(i = 1:dim(path_load)[1], .packages = c("tidyverse"),
                   .combine = rbind) %dopar% {
                     # parameter of the sub-scenarios
                     parameters <- readRDS(paste0(path_load$path_load_data[i], "..parameters.rds"))
                     
                     # loading analyzed data AND selecting only significant trials and where
                     # survival median of control and treatment group can be calculated
                     # for further analysis
                     df.raw.list <- apply(as.matrix(parameters), MARGIN=1, FUN=function(x){
                       rbind(readRDS(paste0(path_load$path_load_ana[i],
                                            "..beta.",     x["beta"],
                                            ".p.C.",       x["p.C"],
                                            ".p.C.C.",     x["p.C.C"],
                                            ".p.C.T.",     x["p.C.T"],
                                            ".r.",         x["r"],
                                            ".med.C.",     x["med.C"],
                                            ".shape.C.T.", x["shape.C.T"],
                                            ".HR.var.",    x["HR.var"],
                                            ".designHR.",  x["designHR"],
                                            ".accrual.",   x["accrual.time"],
                                            ".fu.",        x["fu.time"],
                                            ".rds")) %>%
                               filter(sig %in% c(1,2) & !is.na(median.gain)) %>%
                               mutate(
                                 beta      = x["beta"],
                                 p.C       = x["p.C"],
                                 p.C.C     = x["p.C.C"],
                                 p.C.T     = x["p.C.T"],
                                 r         = x["r"],
                                 med.C     = x["med.C"],
                                 shape.C.T = x["shape.C.T"],
                                 HR.var    = x["HR.var"],
                                 designHR  = x["designHR"],
                                 accrual   = x["accrual.time"],
                                 fu        = x["fu.time"]
                               ) %>%
                               select(beta, p.C, p.C.C, p.C.T, r, med.C, shape.C.T, HR.var,
                                      designHR, accrual, fu,
                                      ESMO, IQWiG_RR, ModIQWiG_HR, ASCO))
                     })
                     # combining needed columns to overall data sets of scenario i
                     df <- bind_rows(df.raw.list) %>%
                       mutate(Scenario = as.character(path_load$Scenario[i]),
                              ESMO.1vs234 = factor(ifelse(ESMO == "1", "1", "234")),
                              ESMO.12vs34 = factor(ifelse(ESMO %in% c("1", "2"), "12", "34")),
                              ESMO.123vs4 = factor(ifelse(ESMO %in% c("1", "2", "3"), "123", "4")),
                              ESMO = factor(ESMO),
                              
                              IQWiG_RR.4vs56 = factor(ifelse(IQWiG_RR == "4", "4", "56")),
                              IQWiG_RR.45vs6 = factor(ifelse(IQWiG_RR %in% c("4", "5"), "45", "6")),
                              IQWiG_RR = factor(IQWiG_RR),
                              
                              ModIQWiG_HR.4vs56 = factor(ifelse(ModIQWiG_HR == "4", "4", "56")),
                              ModIQWiG_HR.45vs6 = factor(ifelse(ModIQWiG_HR %in% c("4", "5"), "45", "6")),
                              ModIQWiG_HR = factor(ModIQWiG_HR)
                       )
                     
                     # Saving results of scenario i
                     if(i==2){ # Scenario 2
                       saveRDS(df %>% filter(HR.var<1),
                               file = paste0(path_saving, "df_cutpoints_", as.character(path_load$Scenario[i]), "Overpowered.rds"))
                       saveRDS(df %>% filter(HR.var>1),
                               file = paste0(path_saving, "df_cutpoints_", as.character(path_load$Scenario[i]), "Underpowered.rds"))
                     }else if(i==3){ # Scenario 3 WEIB
                       saveRDS(df,
                               file = paste0(path_saving, "df_cutpoints_", as.character(path_load$Scenario[i]), ".rds"))
                       saveRDS(df %>% filter(shape.C.T<1),
                               file = paste0(path_saving, "df_cutpoints_", as.character(path_load$Scenario[i]), "DecreasingHazard.rds"))
                       saveRDS(df %>% filter(shape.C.T>1),
                               file = paste0(path_saving, "df_cutpoints_", as.character(path_load$Scenario[i]), "IncreasingHazard.rds"))
                     }else if(i==4){ # Scenario 3 GOMP
                       saveRDS(df,
                               file = paste0(path_saving, "df_cutpoints_", as.character(path_load$Scenario[i]), ".rds"))
                       saveRDS(df %>% filter(shape.C.T<0),
                               file = paste0(path_saving, "df_cutpoints_", as.character(path_load$Scenario[i]), "DecreasingHazard.rds"))
                       saveRDS(df %>% filter(shape.C.T>0),
                               file = paste0(path_saving, "df_cutpoints_", as.character(path_load$Scenario[i]), "IncreasingHazard.rds"))
                     }else if(i==6){# Scenario 5
                       saveRDS(df %>% filter(r<1),
                               file = paste0(path_saving, "df_cutpoints_", as.character(path_load$Scenario[i]), "RSmaller1.rds"))
                       saveRDS(df %>% filter(r>1),
                               file = paste0(path_saving, "df_cutpoints_", as.character(path_load$Scenario[i]), "RGreater1.rds"))
                     }else if(i==8){# Scenario 7
                       saveRDS(df %>% filter(p.C.C==0.2),
                               file = paste0(path_saving, "df_cutpoints_", as.character(path_load$Scenario[i]), "CtlSmaller.rds"))
                       saveRDS(df %>% filter(p.C.C==0.4),
                               file = paste0(path_saving, "df_cutpoints_", as.character(path_load$Scenario[i]), "CtlGreater.rds"))
                     }else{
                       saveRDS(df,
                               file = paste0(path_saving, "df_cutpoints_", as.character(path_load$Scenario[i]), ".rds"))
                     }
                     rm(df)
                     rm(df.raw.list)
                   }
# Stopping parallel computing 
stopCluster(cl)

(time.1 <- Sys.time() - start.time.1)


#--------------------- Analyzing combined data sets ----------------------------
################################################################################
######################### Cutoff Values ########################################
################################################################################

### functions
SvenssonCutoff <- function(df., numericVariable = "ASCO", categoricalVariable){
  cutoff <- rep(NA, length(unique(df.[,categoricalVariable]))-1)
  for (i in 1:(length(unique(df.[,categoricalVariable]))-1)) {
    if(i==1){
      cutoff[i] <-
        uniroot(f = function(x){
          table(df.[,categoricalVariable], exclude = levels(df.[,categoricalVariable])[!levels(df.[,categoricalVariable]) %in% as.character(unique(df.[,categoricalVariable]))])[[i]] - sum(df.[,numericVariable]<=x)
        },
        interval = c(0,100000),tol = .Machine$double.eps^0.35)$root
    }else{
      cutoff[i] <-
        uniroot(f = function(x){
          table(df.[,categoricalVariable], exclude = levels(df.[,categoricalVariable])[!levels(df.[,categoricalVariable]) %in% as.character(unique(df.[,categoricalVariable]))])[[i]] - sum(df.[,numericVariable]<=x & df.[,numericVariable]>cutoff[i-1])
        },
        interval = c(cutoff[i-1],100000), tol = .Machine$double.eps^0.35)$root
    }
  }
  return(cutoff)
}

WeightedCohensKappa_ESMO_ASCO  <- function(x, df.){
  # df. has to be a data.table (so that the data steps are faster)!
  
  # order cutoff values, so that it stays logical
  cutoffs <- sort(c(x[1], x[2], x[3]))
  
  # Compute ASCO categories using the cutoff values from "cutoffs"
  df.[, ASCO.catESMO := factor(
    ifelse(ASCO<=cutoffs[1], "1",
           ifelse(ASCO>cutoffs[1] & ASCO<=cutoffs[2], "2",
                  ifelse(ASCO>cutoffs[2] & ASCO<=cutoffs[3], "3", "4"))),
    levels = c(1,2,3,4))]
  
  # Compute table
  ASCO_ESMO_table <- df.[, table(ESMO, ASCO.catESMO, exclude=NULL)]
  
  # Calculating weighted Cohen Kappa (with quadratic weights)
  WeightedCohenKappa <- Kappa(ASCO_ESMO_table, weights = "Fleiss-Cohen")$Weighted[[1]]
  
  # Returning 99999 if WeightedCohenKappa is negative, otherwise
  # returning negative WeightedCohenKappa (because the optimizer function "optim"
  # performs minimization)
  return(ifelse(WeightedCohenKappa<0, 99999, -WeightedCohenKappa))
}
WeightedCohensKappa_IQWiG_RR_ASCO  <- function(x, df.){
  # df. has to be a data.table (so that the data steps are faster)!
  
  # order cutoff values, so that it stays logical
  cutoffs <- sort(c(x[1], x[2]))
  
  # Compute ASCO categories using the cutoff values from "cutoffs"
  df.[, ASCO.catIQWiG := factor(
    ifelse(ASCO<=cutoffs[1], "4",
           ifelse(ASCO>cutoffs[1] & ASCO<=cutoffs[2], "5", "6")),
    levels = c(4,5,6))]
  
  # Compute table
  ASCO_IQWiG_table <- df.[, table(IQWiG_RR, ASCO.catIQWiG, exclude=NULL)]
  
  # Calculating weighted Cohen Kappa (with quadratic weights)
  WeightedCohenKappa <- Kappa(ASCO_IQWiG_table, weights = "Fleiss-Cohen")$Weighted[[1]]
  
  # Returning 99999 if WeightedCohenKappa is negative, otherwise
  # returning negative WeightedCohenKappa (because the optimizer function "optim"
  # performs minimization)
  return(ifelse(WeightedCohenKappa<0, 99999, -WeightedCohenKappa))
}
WeightedCohensKappa_ModIQWiG_HR_ASCO  <- function(x, df.){
  # df. has to be a data.table (so that the data steps are faster)!
  
  # order cutoff values, so that it stays logical
  cutoffs <- sort(c(x[1], x[2]))
  
  # Compute ASCO categories using the cutoff values from "cutoffs"
  df.[, ASCO.catModIQWiG_HR := factor(
    ifelse(ASCO<=cutoffs[1], "4",
           ifelse(ASCO>cutoffs[1] & ASCO<=cutoffs[2], "5", "6")),
    levels = c(4,5,6))]
  
  # Compute table
  ASCO_ModIQWiG_HR_table <- df.[, table(ModIQWiG_HR, ASCO.catModIQWiG_HR, exclude=NULL)]
  
  # Calculating weighted Cohen Kappa (with quadratic weights)
  WeightedCohenKappa <- Kappa(ASCO_ModIQWiG_HR_table, weights = "Fleiss-Cohen")$Weighted[[1]]
  
  # Returning 99999 if WeightedCohenKappa is negative, otherwise
  # returning negative WeightedCohenKappa (because the optimizer function "optim"
  # performs minimization)
  return(ifelse(WeightedCohenKappa<0, 99999, -WeightedCohenKappa))
}


# Calculating cutpoints
results_cutpoints <- matrix(
  list(), nrow=15, ncol = 21,
  dimnames = list(Scenarios = c("Scen1", "Scen2Overpowered", "Scen2Underpowered",
                                "Scen3WEIB", "Scen3WEIBDecreasingHazard", "Scen3WEIBIncreasingHazard",
                                "Scen3GOMP", "Scen3GOMPDecreasingHazard", "Scen3GOMPIncreasingHazard",
                                "Scen4",
                                "Scen5RSmaller1", "Scen5RGreater1", "Scen6",
                                "Scen7CtlSmaller", "Scen7CtlGreater"),
                  MethodsForCutOffs = c(
                    "ESMO_CohensKappa1","ESMO_CohensKappa12","ESMO_CohensKappa123",
                    "ESMO_roc01_1", "ESMO_roc01_12", "ESMO_roc01_123",
                    "ESMO_Svensson1", "ESMO_Svensson12", "ESMO_Svensson123",

                    "IQWiG_RR_CohensKappa4","IQWiG_RR_CohensKappa45",
                    "IQWiG_RR_roc01_4", "IQWiG_RR_roc01_45",
                    "IQWiG_RR_Svensson4", "IQWiG_RR_Svensson45",
                    
                    "ModIQWiG_HR_CohensKappa4","ModIQWiG_HR_CohensKappa45",
                    "ModIQWiG_HR_roc01_4", "ModIQWiG_HR_roc01_45",
                    "ModIQWiG_HR_Svensson4", "ModIQWiG_HR_Svensson45")))

for (i in rownames(results_cutpoints)) {
  ####### Loading data
  df <- readRDS(paste0(path_saving, "df_cutpoints_", i, ".rds"))
  
  ####### Calculating ASCO cutoff values for ESMO
  #### ROC: roc01
  results_cutpoints[i, "ESMO_roc01_1"][[1]] <-
    cutpointr(data = df, x = ASCO, class = ESMO.1vs234,
              pos_class = "234", neg_class = "1", direction = ">=",
              method = minimize_metric, metric = roc01, break_ties = c) %>%
    select(optimal_cutpoint, roc01) %>%
    unnest(cols = c(optimal_cutpoint, roc01))
  
  results_cutpoints[i, "ESMO_roc01_12"][[1]] <-
    cutpointr(data = df,x = ASCO, class = ESMO.12vs34,
              pos_class = "34", neg_class = "12", direction = ">=",
              method = minimize_metric, metric = roc01, break_ties = c) %>%
    select(optimal_cutpoint, roc01) %>%
    unnest(cols = c(optimal_cutpoint, roc01))
  
  results_cutpoints[i, "ESMO_roc01_123"][[1]] <-
    cutpointr(data = df,x = ASCO, class = ESMO.123vs4,
              pos_class = "4", neg_class = "123", direction = ">=",
              method = minimize_metric, metric = roc01, break_ties = c) %>%
    select(optimal_cutpoint, roc01) %>%
    unnest(cols = c(optimal_cutpoint, roc01))
  
  #### Svenssons method
  SvenssonsCutoff_ESMO <- SvenssonCutoff(df. = as.data.frame(df), numericVariable = "ASCO", categoricalVariable = "ESMO")
  results_cutpoints[i, "ESMO_Svensson1"][[1]]   <- SvenssonsCutoff_ESMO[1]
  results_cutpoints[i, "ESMO_Svensson12"][[1]]  <- SvenssonsCutoff_ESMO[2]
  results_cutpoints[i, "ESMO_Svensson123"][[1]] <- SvenssonsCutoff_ESMO[3]
  
  #### Maximizing cohens kappa approach
  # ## Checking convexity (if convex starting point would not matter)
  # set.seed(123456789)
  # x1 <- runif(n=100, min = min(df$ASCO), max = max(df$ASCO))
  # x2 <- runif(n=100, min = min(df$ASCO), max = max(df$ASCO))
  # x3 <- runif(n=100, min = min(df$ASCO), max = max(df$ASCO))
  #
  # y1 <- runif(n=100, min = min(df$ASCO), max = max(df$ASCO))
  # y2 <- runif(n=100, min = min(df$ASCO), max = max(df$ASCO))
  # y3 <- runif(n=100, min = min(df$ASCO), max = max(df$ASCO))
  #
  # theta <- seq(0,1,by = 0.2)
  #
  # res_j <- sapply(theta, function(theta.i, x, y){
  #    # Sorting cutoff values 
  #    x <- sort(x)
  #    y <- sort(y)
  #
  #    # Calculating left and right side
  #    LeftSide <- WeightedCohensKappa_ESMO_ASCO(df. = setDT(df)[, list(ASCO, ESMO)], x = c(theta.i*x[1] + (1-theta.i)*y[1],
  #                                                                                               theta.i*x[2] + (1-theta.i)*y[2],
  #                                                                                               theta.i*x[3] + (1-theta.i)*y[3]))
  #    RightSide <- theta.i*WeightedCohensKappa_ESMO_ASCO(df. = setDT(df)[, list(ASCO, ESMO)], x = c(x[1], x[2], x[3])) +
  #       (1-theta.i)*WeightedCohensKappa_ESMO_ASCO(df. = setDT(df)[, list(ASCO, ESMO)], x = c(y[1], y[2], y[3]))
  #    # Checking convexity
  #    return(LeftSide<=RightSide)
  # }, x = c(x1[1], x2[1], x3[1]), y = c(y1[1], y2[1], y3[1]))
  #
  #
  # # Example that of no convexity
  # theta.i <- 0.2
  # x <- c(x1[1], x2[1], x3[1])[order(c(x1[1], x2[1], x3[1]))]
  # y <- c(y1[1], y2[1], y3[1])[order(c(y1[1], y2[1], y3[1]))]
  # LeftSide <- WeightedCohensKappa_ESMO_ASCO(df. = setDT(df)[, list(ASCO, ESMO)], x = c(theta.i*x[1] + (1-theta.i)*y[1],
  #                                                                                            theta.i*x[2] + (1-theta.i)*y[2],
  #                                                                                            theta.i*x[3] + (1-theta.i)*y[3]))
  # RightSide <- theta.i*WeightedCohensKappa_ESMO_ASCO(df. = setDT(df)[, list(ASCO, ESMO)], x = c(x[1], x[2], x[3])) +
  #    (1-theta.i)*WeightedCohensKappa_ESMO_ASCO(df. = setDT(df)[, list(ASCO, ESMO)], x = c(y[1], y[2], y[3]))
  # LeftSide<=RightSide
  # LeftSide
  # RightSide
  # #--> no convexity --> multiple minima possible --> using various starting points
  
  
  ## Local Minima with different start values (--> parallel computing)
  n_start_values <- 40
  set.seed(123456)
  start_values <- data.frame(matrix(NA, nrow = n_start_values, ncol = 4))
  start_values[,1] <- c(
    c("roc01+Paper", "Svensson+Paper", "roc01", "Svensson", "MeanOfGroups"),
    rep("Random", n_start_values-5))
  start_values[1:5,2:4]  <- matrix(c(
    # roc01 + Paper
    mean(results_cutpoints[i, "ESMO_roc01_1"][[1]]$optimal_cutpoint), 42.250, 44.815,
    # Svensson + Paper
    results_cutpoints[i, "ESMO_Svensson1"][[1]], 42.250, 44.815,
    # roc01
    mean(results_cutpoints[i, "ESMO_roc01_1"][[1]]$optimal_cutpoint),
    mean(results_cutpoints[i, "ESMO_roc01_12"][[1]]$optimal_cutpoint),
    mean(results_cutpoints[i, "ESMO_roc01_123"][[1]]$optimal_cutpoint),
    # Svensson
    results_cutpoints[i, "ESMO_Svensson1"][[1]],
    results_cutpoints[i, "ESMO_Svensson12"][[1]],
    results_cutpoints[i, "ESMO_Svensson123"][[1]],
    # mean of group means
    mean(c(mean(df$ASCO[which(df$ESMO=="1")]), mean(df$ASCO[which(df$ESMO=="2")]))),
    mean(c(mean(df$ASCO[which(df$ESMO=="2")]), mean(df$ASCO[which(df$ESMO=="3")]))),
    mean(c(mean(df$ASCO[which(df$ESMO=="3")]), mean(df$ASCO[which(df$ESMO=="4")])))
  ),
  ncol = 3, byrow = TRUE)
  
  for (j in 6:n_start_values) {
    start_values[j,2:4] <- sort(runif(n = 3, min = min(df$ASCO), max = max(df$ASCO)))
  }
  
  # Saving data as "data.table", selecting only needed coloumns and keep df itsel as "data.frame"
  df_help <- setDT(df)[, list(ASCO, ESMO)]
  df <- as.data.frame(df)
  
  start.time.1 <- Sys.time()
  # Setting up parallel computing
  cl <- makeCluster(num.cl) 
  registerDoParallel(cl)
  
  # Starting parallel computing
  MaxKappa_ESMO_ASCO <- foreach(
    j = iter(start_values, by='row'), .packages = c("vcd", "data.table"), 
    .combine = rbind) %dopar% {
      opt_NelderMead <- optim(fn=WeightedCohensKappa_ESMO_ASCO, par = j[2:4],
                              df. = df_help,
                              method = "Nelder-Mead",
                              control = list(
                                # Relative convergence tolerance:
                                # The algorithm stops if it is unable to reduce the value
                                # by a factor of reltol * (abs(val) + reltol) at a step.
                                # Defaults to sqrt(.Machine$double.eps), typically about 1e-8.
                                reltol = sqrt(.Machine$double.eps)#,
                                #trace = 6
                              ))
      
      return <- data.frame(start_value_name = j[1],
                           start_value1 = j[2],
                           start_value2 = j[3],
                           start_value3 = j[4],
                           optimal_value1 = sort(opt_NelderMead$par)[1],
                           optimal_value2 = sort(opt_NelderMead$par)[2],
                           optimal_value3 = sort(opt_NelderMead$par)[3],
                           Converge       = opt_NelderMead$convergence,
                           KappaValue     = opt_NelderMead$value
      )
      colnames(return) <- c("start_values_name", "start_values1",
                            "start_values2", "start_values3",
                            colnames(return)[5:9])
      return(return)
    }
  # Stopping parallel computing
  stopCluster(cl)
  
  (time.1 <- Sys.time() - start.time.1)
  
  # Order maximizing Cohens Kappa results / output
  MaxKappa_ESMO_ASCO_orderd <- MaxKappa_ESMO_ASCO[order(MaxKappa_ESMO_ASCO$KappaValue),]
  
  # get optimal cutoffs and put them into "results_cutpoints"
  results_cutpoints[i, "ESMO_CohensKappa1"][[1]]   <- MaxKappa_ESMO_ASCO_orderd[1,c(5,9)]
  results_cutpoints[i, "ESMO_CohensKappa12"][[1]]  <- MaxKappa_ESMO_ASCO_orderd[1,c(6,9)]
  results_cutpoints[i, "ESMO_CohensKappa123"][[1]] <- MaxKappa_ESMO_ASCO_orderd[1,c(7,9)]
  
  # Saving maximizing Cohen Kappa results (severally because that takes the longest)
  saveRDS(MaxKappa_ESMO_ASCO_orderd,
          file = paste0(path_saving, "results_MaxKappa_ESMO_ASCO_orderd_", i, ".rds"))
  
  
  ####### Calculating ASCO cutoff values for IQWiG_RR
  #### ROC: roc01
  results_cutpoints[i, "IQWiG_RR_roc01_4"][[1]] <-
    cutpointr(data = df, x = ASCO, class = IQWiG_RR.4vs56,
              pos_class = "56", neg_class = "4", direction = ">=",
              method = minimize_metric, metric = roc01, break_ties = c) %>%
    select(optimal_cutpoint, roc01) %>%
    unnest(cols = c(optimal_cutpoint, roc01))
  
  results_cutpoints[i, "IQWiG_RR_roc01_45"][[1]] <-
    cutpointr(data = df,x = ASCO, class = IQWiG_RR.45vs6,
              pos_class = "6", neg_class = "45", direction = ">=",
              method = minimize_metric, metric = roc01, break_ties = c) %>%
    select(optimal_cutpoint, roc01) %>%
    unnest(cols = c(optimal_cutpoint, roc01))
  
  #### Svenssons method
  SvenssonsCutoff_IQWiG_RR <- SvenssonCutoff(df. = as.data.frame(df), numericVariable = "ASCO", categoricalVariable = "IQWiG_RR")
  results_cutpoints[i, "IQWiG_RR_Svensson4"][[1]]   <- SvenssonsCutoff_IQWiG_RR[1]
  results_cutpoints[i, "IQWiG_RR_Svensson45"][[1]]  <- SvenssonsCutoff_IQWiG_RR[2]
  
  
  #### Maximizing cohens kappa approach
  # ## Checking convexity (if convex starting point would not matter)
  # set.seed(123456789)
  # x1 <- runif(n=100, min = min(df$ASCO), max = max(df$ASCO))
  # x2 <- runif(n=100, min = min(df$ASCO), max = max(df$ASCO))
  #
  # y1 <- runif(n=100, min = min(df$ASCO), max = max(df$ASCO))
  # y2 <- runif(n=100, min = min(df$ASCO), max = max(df$ASCO))
  #
  # theta <- seq(0,1,by = 0.2)
  #
  # res_j <- sapply(theta, function(theta.i, x, y){
  #    # Sorting cutoff values 
  #    x <- sort(x)
  #    y <- sort(y)
  #
  #    # Calculating left and right side
  #    LeftSide <- WeightedCohensKappa_IQWiG_RR_ASCO(df. = setDT(df)[, list(ASCO, IQWiG_RR)], x = c(theta.i*x[1] + (1-theta.i)*y[1],
  #                                                                                                 theta.i*x[2] + (1-theta.i)*y[2]))
  #    RightSide <- theta.i*WeightedCohensKappa_IQWiG_RR_ASCO(df. = setDT(df)[, list(ASCO, IQWiG_RR)], x = c(x[1], x[2])) +
  #       (1-theta.i)*WeightedCohensKappa_IQWiG_ASCO(df. = setDT(df)[, list(ASCO, IQWiG_RR)], x = c(y[1], y[2]))
  #    # Checking convexity
  #    return(LeftSide<=RightSide)
  # }, x = c(x1[1], x2[1]), y = c(y1[1], y2[1]))
  #
  #
  # # Example that of no convexity
  # theta.i <- 0.2
  # x <- c(x1[1], x2[1])[order(c(x1[1], x2[1]))]
  # y <- c(y1[1], y2[1])[order(c(y1[1], y2[1]))]
  # LeftSide <- WeightedCohensKappa_IQWiG_RR_ASCO(df. = setDT(df)[, list(ASCO, IQWiG_RR)], x = c(theta.i*x[1] + (1-theta.i)*y[1],
  #                                                                                            theta.i*x[2] + (1-theta.i)*y[2]))
  # RightSide <- theta.i*WeightedCohensKappa_IQWiG_RR_ASCO(df. = setDT(df)[, list(ASCO, IQWiG_RR)], x = c(x[1], x[2])) +
  #    (1-theta.i)*WeightedCohensKappa_IQWiG_RR_ASCO(df. = setDT(df)[, list(ASCO, IQWiG_RR)], x = c(y[1], y[2]))
  # LeftSide<=RightSide
  # LeftSide
  # RightSide
  # #--> no convexity --> multiple minima possible --> using various starting points
  
  ## Local Minima with different start values (--> parallel computing)
  n_start_values <- 40
  set.seed(123456)
  start_values <- data.frame(matrix(NA, nrow = n_start_values, ncol = 3))
  start_values[,1] <- c(
    c("roc01", "Svensson", "MeanOfGroups"),
    rep("Random", n_start_values-3))
  start_values[1:3,2:3]  <- matrix(c(
    # roc01
    mean(results_cutpoints[i, "IQWiG_RR_roc01_4"][[1]]$optimal_cutpoint), mean(results_cutpoints[i, "IQWiG_RR_roc01_45"][[1]]$optimal_cutpoint),
    # svenson
    results_cutpoints[i, "IQWiG_RR_Svensson4"][[1]], results_cutpoints[i, "IQWiG_RR_Svensson45"][[1]],
    # mean of group means
    mean(c(mean(df$ASCO[which(df$IQWiG_RR=="4")]), mean(df$ASCO[which(df$IQWiG_RR=="5")]))),
    mean(c(mean(df$ASCO[which(df$IQWiG_RR=="5")]), mean(df$ASCO[which(df$IQWiG_RR=="6")])))
  ),
  ncol = 2, byrow = TRUE)
  
  for (j in 4:n_start_values) {
    start_values[j,2:3] <- sort(runif(n = 2, min = min(df$ASCO), max = max(df$ASCO)))
  }
  
  # Saving data as "data.table", selecting only needed coloumns and keep df itsel as "data.frame"
  df_help <- setDT(df)[, list(ASCO, IQWiG_RR)]
  df <- as.data.frame(df)
  
  start.time.1 <- Sys.time()
  # Setting up parallel computing
  cl <- makeCluster(num.cl) 
  registerDoParallel(cl)
  
  # starting parallel computing
  MaxKappa_IQWiG_RR_ASCO <- foreach(
    j = iter(start_values, by='row'), .packages = c("vcd", "data.table"),
    .combine = rbind) %dopar% {
      opt_NelderMead <- optim(fn=WeightedCohensKappa_IQWiG_RR_ASCO, par = j[2:3],
                              df. = df_help,
                              method = "Nelder-Mead",
                              control = list(
                                # Relative convergence tolerance:
                                # The algorithm stops if it is unable to reduce the value
                                # by a factor of reltol * (abs(val) + reltol) at a step.
                                # Defaults to sqrt(.Machine$double.eps), typically about 1e-8.
                                reltol = sqrt(.Machine$double.eps)#,
                                #trace = 6
                              ))
      
      return <- data.frame(start_values_name = j[1],
                           start_value1 = j[2],
                           start_value2 = j[3],
                           optimal_value1 = sort(opt_NelderMead$par)[1],
                           optimal_value2 = sort(opt_NelderMead$par)[2],
                           Converge       = opt_NelderMead$convergence,
                           KappaValue     = opt_NelderMead$value)
      colnames(return) <- c("start_values_name", "start_values1", "start_values2",
                            colnames(return)[4:7])
      return(return)
    }
  # Stopping parallel computing
  stopCluster(cl)
  
  (time.1 <- Sys.time() - start.time.1)
  
  # Order maximizing Cohens Kappa results / output
  MaxKappa_IQWiG_RR_ASCO_orderd <- MaxKappa_IQWiG_RR_ASCO[order(MaxKappa_IQWiG_RR_ASCO$KappaValue),]
  
  # get optimal cutoffs and put them into "results_cutpoints"
  results_cutpoints[i, "IQWiG_RR_CohensKappa4"][[1]]  <- MaxKappa_IQWiG_RR_ASCO_orderd[1, c(4,7)]
  results_cutpoints[i, "IQWiG_RR_CohensKappa45"][[1]] <- MaxKappa_IQWiG_RR_ASCO_orderd[1, c(5,7)]
  
  # Saving maximizing Cohen Kappa results (severally because that takes the longest)
  saveRDS(MaxKappa_IQWiG_RR_ASCO_orderd,
          file = paste0(path_saving, "results_MaxKappa_IQWiG_RR_ASCO_orderd_", i, ".rds"))
  
  
  
  
  ####### Calculating ASCO cutoff values for ModIQWiG_HR
  #### ROC: roc01
  results_cutpoints[i, "ModIQWiG_HR_roc01_4"][[1]] <-
    cutpointr(data = df, x = ASCO, class = ModIQWiG_HR.4vs56,
              pos_class = "56", neg_class = "4", direction = ">=",
              method = minimize_metric, metric = roc01, break_ties = c) %>%
    select(optimal_cutpoint, roc01) %>%
    unnest(cols = c(optimal_cutpoint, roc01))
  
  results_cutpoints[i, "ModIQWiG_HR_roc01_45"][[1]] <-
    cutpointr(data = df,x = ASCO, class = ModIQWiG_HR.45vs6,
              pos_class = "6", neg_class = "45", direction = ">=",
              method = minimize_metric, metric = roc01, break_ties = c) %>%
    select(optimal_cutpoint, roc01) %>%
    unnest(cols = c(optimal_cutpoint, roc01))
  
  #### Svenssons method
  SvenssonsCutoff_ModIQWiG_HR <- SvenssonCutoff(df. = as.data.frame(df), numericVariable = "ASCO", categoricalVariable = "ModIQWiG_HR")
  results_cutpoints[i, "ModIQWiG_HR_Svensson4"][[1]]   <- SvenssonsCutoff_ModIQWiG_HR[1]
  results_cutpoints[i, "ModIQWiG_HR_Svensson45"][[1]]  <- SvenssonsCutoff_ModIQWiG_HR[2]
  
  
  #### Maximizing cohens kappa approach
  # ## Checking convexity (if convex starting point would not matter)
  # set.seed(123456789)
  # x1 <- runif(n=100, min = min(df$ASCO), max = max(df$ASCO))
  # x2 <- runif(n=100, min = min(df$ASCO), max = max(df$ASCO))
  #
  # y1 <- runif(n=100, min = min(df$ASCO), max = max(df$ASCO))
  # y2 <- runif(n=100, min = min(df$ASCO), max = max(df$ASCO))
  #
  # theta <- seq(0,1,by = 0.2)
  #
  # res_j <- sapply(theta, function(theta.i, x, y){
  #    # Sorting cutoff values 
  #    x <- sort(x)
  #    y <- sort(y)
  #
  #    # Calculating left and right side
  #    LeftSide <- WeightedCohensKappa_ModIQWiG_HR_ASCO(df. = setDT(df)[, list(ASCO, ModIQWiG_HR)], x = c(theta.i*x[1] + (1-theta.i)*y[1],
  #                                                                                                 theta.i*x[2] + (1-theta.i)*y[2]))
  #    RightSide <- theta.i*WeightedCohensKappa_ModIQWiG_HR_ASCO(df. = setDT(df)[, list(ASCO, ModIQWiG_HR)], x = c(x[1], x[2])) +
  #       (1-theta.i)*WeightedCohensKappa_ModIQWiG_HR_ASCO(df. = setDT(df)[, list(ASCO, ModIQWiG_HR)], x = c(y[1], y[2]))
  #    # Checking convexity
  #    return(LeftSide<=RightSide)
  # }, x = c(x1[1], x2[1]), y = c(y1[1], y2[1]))
  #
  #
  # # Example that of no convexity
  # theta.i <- 0.2
  # x <- c(x1[1], x2[1])[order(c(x1[1], x2[1]))]
  # y <- c(y1[1], y2[1])[order(c(y1[1], y2[1]))]
  # LeftSide <- WeightedCohensKappa_ModIQWiG_HR_ASCO(df. = setDT(df)[, list(ASCO, ModIQWiG_HR)], x = c(theta.i*x[1] + (1-theta.i)*y[1],
  #                                                                                            theta.i*x[2] + (1-theta.i)*y[2]))
  # RightSide <- theta.i*WeightedCohensKappa_ModIQWiG_HR_ASCO(df. = setDT(df)[, list(ASCO, ModIQWiG_HR)], x = c(x[1], x[2])) +
  #    (1-theta.i)*WeightedCohensKappa_ModIQWiG_HR_ASCO(df. = setDT(df)[, list(ASCO, ModIQWiG_HR)], x = c(y[1], y[2]))
  # LeftSide<=RightSide
  # LeftSide
  # RightSide
  # #--> no convexity --> multiple minima possible --> using various starting points
  
  ## Local Minima with different start values (--> parallel computing)
  n_start_values <- 40
  set.seed(123456)
  start_values <- data.frame(matrix(NA, nrow = n_start_values, ncol = 3))
  start_values[,1] <- c(
    c("roc01", "Svensson", "MeanOfGroups"),
    rep("Random", n_start_values-3))
  start_values[1:3,2:3]  <- matrix(c(
    # roc01
    mean(results_cutpoints[i, "ModIQWiG_HR_roc01_4"][[1]]$optimal_cutpoint), mean(results_cutpoints[i, "ModIQWiG_HR_roc01_45"][[1]]$optimal_cutpoint),
    # Svensson
    results_cutpoints[i, "ModIQWiG_HR_Svensson4"][[1]], results_cutpoints[i, "ModIQWiG_HR_Svensson45"][[1]],
    # mean of group means
    mean(c(mean(df$ASCO[which(df$ModIQWiG_HR=="4")]), mean(df$ASCO[which(df$ModIQWiG_HR=="5")]))),
    mean(c(mean(df$ASCO[which(df$ModIQWiG_HR=="5")]), mean(df$ASCO[which(df$ModIQWiG_HR=="6")])))
  ),
  ncol = 2, byrow = TRUE)
  
  for (j in 4:n_start_values) {
    start_values[j,2:3] <- sort(runif(n = 2, min = min(df$ASCO), max = max(df$ASCO)))
  }
  
  # Saving data as "data.table", selecting only needed coloumns and keep df itsel as "data.frame"
  df_help <- setDT(df)[, list(ASCO, ModIQWiG_HR)]
  df <- as.data.frame(df)
  
  start.time.1 <- Sys.time()
  # Setting up parallel computing
  cl <- makeCluster(num.cl) 
  registerDoParallel(cl)
  
  # starting parallel computing
  MaxKappa_ModIQWiG_HR_ASCO <- foreach(
    j = iter(start_values, by='row'), .packages = c("vcd", "data.table"),
    .combine = rbind) %dopar% {
      opt_NelderMead <- optim(fn=WeightedCohensKappa_ModIQWiG_HR_ASCO, par = j[2:3],
                              df. = df_help,
                              method = "Nelder-Mead",
                              control = list(
                                # Relative convergence tolerance:
                                # The algorithm stops if it is unable to reduce the value
                                # by a factor of reltol * (abs(val) + reltol) at a step.
                                # Defaults to sqrt(.Machine$double.eps), typically about 1e-8.
                                reltol = sqrt(.Machine$double.eps)#,
                                #trace = 6
                              ))
      return <- data.frame(start_values_name = j[1],
                           start_value1 = j[2],
                           start_value2 = j[3],
                           optimal_value1 = sort(opt_NelderMead$par)[1],
                           optimal_value2 = sort(opt_NelderMead$par)[2],
                           Converge       = opt_NelderMead$convergence,
                           KappaValue     = opt_NelderMead$value)
      colnames(return) <- c("start_values_name", "start_values1", "start_values2",
                            colnames(return)[4:7])
      return(return)
    }
  # Stopping parallel computing
  stopCluster(cl)
  
  (time.1 <- Sys.time() - start.time.1)
  
  # Order maximizing Cohens Kappa results / output
  MaxKappa_ModIQWiG_HR_ASCO_orderd <- MaxKappa_ModIQWiG_HR_ASCO[order(MaxKappa_ModIQWiG_HR_ASCO$KappaValue),]
  
  # get optimal cutoffs and put them into "results_cutpoints"
  results_cutpoints[i, "ModIQWiG_HR_CohensKappa4"][[1]]  <- MaxKappa_ModIQWiG_HR_ASCO_orderd[1, c(4,7)]
  results_cutpoints[i, "ModIQWiG_HR_CohensKappa45"][[1]] <- MaxKappa_ModIQWiG_HR_ASCO_orderd[1, c(5,7)]
  
  # Saving maximizing Cohen Kappa results (severally because that takes the longest)
  saveRDS(MaxKappa_ModIQWiG_HR_ASCO_orderd,
          file = paste0(path_saving, "results_MaxKappa_ModIQWiG_HR_ASCO_orderd_", i, ".rds"))
}

### saving the complete results
saveRDS(results_cutpoints,
        file = paste0(path_saving, "results_cutpoints.rds"))