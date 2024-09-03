################################################################################
#########       Analysis of the generated data of Simulation 2         #########
################################################################################ 

# Setting working directory (where functions and Simulation.seed are saved)
setwd("") 

#---------------------- Needed Packages and Functions --------------------------
library(tidyverse)
library(foreach)
library(doParallel)

source("Sim2_0_Functions_Analysis.R")

#-------------------- Same parameter in all scenarios  -------------------------
# Number of used cores for analysis
num.cl <- 40

# Loading / Saving paths
path_data <- "" # e.g. ".../Simulation2/Data/" 
path_ana  <- "" # e.g. ".../Simulation2/Ana/"

#------------------------- Standard Scenario 1 ---------------------------------
func.DataAnalysis.Scenario(
  path_load = paste0(path_data, "Scen1"),
  path_save = paste0(path_ana, "Scen1"),
  Scenario = "Scen1"
)

#------------------------------ Scenario 2 -------------------------------------
#---- Influence of incorrectly assumed designHR for sample size calculation ----
func.DataAnalysis.Scenario(
  path_load = paste0(path_data, "Scen2"),
  path_save = paste0(path_ana, "Scen2"),
  Scenario = "Scen2"
)
#------------------------------ Scenario 3 -------------------------------------
#--- Influence of different underlying failure time distributions (WEIBULL) ----
func.DataAnalysis.Scenario(
  path_load = paste0(path_data, "Scen3WEIB"),
  path_save = paste0(path_ana, "Scen3WEIB"),
  Scenario = "Scen3WEIB"
)
#------------------------------ Scenario 3 -------------------------------------
#--- Influence of different underlying failure time distributions (GOMPERTZ) ---
func.DataAnalysis.Scenario(
  path_load = paste0(path_data, "Scen3GOMP"),
  path_save = paste0(path_ana, "Scen3GOMP"),
  Scenario = "Scen3GOMP"
)
#------------------------------ Scenario 4 -------------------------------------
#--- Influence of non-proportional hazards using late treatment effects for ----
#-------------------------- the treatment group --------------------------------
func.DataAnalysis.Scenario(
  path_load = paste0(path_data, "Scen4"),
  path_save = paste0(path_ana, "Scen4"),
  Scenario = "Scen4"
)

#---------------- Combining everything into one final data set -----------------
CombinedResults <- rbind(
  readRDS(paste0(path_ana, "Scen1..CombinedResults.rds")),
  readRDS(paste0(path_ana, "Scen2..CombinedResults.rds")),
  readRDS(paste0(path_ana, "Scen3WEIB..CombinedResults.rds")),
  readRDS(paste0(path_ana, "Scen3GOMP..CombinedResults.rds")),
  readRDS(paste0(path_ana, "Scen4..CombinedResults.rds"))
)
saveRDS(CombinedResults,
        file = paste0(path_ana, "CombinedResults.rds")
)

