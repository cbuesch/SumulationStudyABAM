# SimulationStudyABAM
This github respository provides the R-code to reproduce the results of the PhD-Thesis "Comparison of different additional benefit assessment methods for oncology treatments" of Christopher Alexander Büsch. The R-Code of both simulation studies are shared in the folder "RPrograms". Please read the R-Code instruction below before performing it!

## R-Code: 
The following R-Packages need to be installed before R-Code usage: tidyverse, stringr, cutpointr, vcd, data.table, survival, flexsurv, irr, foreach and doParallel. This can be done by using the code provided hereafter:
```r
install.packages(c("tidyverse", "stringr", "cutpointr", "vcd", "data.table", "survival", "flexsurv", "irr", "foreach", "doParallel"))
```
The provided programs are using the doParallel package for parallel computing in windows systems. In case a unix-like system is used, the package doMC instead of doParallel needs to be installed. Furthermore, all occurrences of 
```r
cl <- makeCluster(num.cl)
registerDoParallel(cl)
```
need to be replaced by 
```r
registerDoMC(num.cl)
```
Moreover, all lines containing 
```r
stopCluster(cl)
```
need to be removed. In addition, the running time of the programs (especially the data generation and data analysis) can take even with multiple cores for parallel computing several days. Therefore, to reduce the running time, 40 cores or more are recommended. Another solution is to reduce the number of iterations (n.sim), which decreases, however, the precision of the simulation study.
In the following, several R-Code-Scripts are explained, which can be found in the folder "RPrograms":
- Costume functions for data generation and analysis of both simulation studies ("Sim1_0_Functions_DataGeneration.R", "Sim2_0_Functions_DataGeneration.R", "Sim1_0_Functions_Analysis.R" and "Sim2_0_Functions_Analysis.R"):
  These four scripts contain all needed functions for the data generation and analysis. Hence, they need to be loaded before conducting the simulations. 
  
- Data generation of both simulation studies ("Sim1_1_DataGeneration.R" and "Sim2_1_DataGeneration.R"):
  These scripts conduct the data generation for all scenarios of both simulation studies. At the beginning the working directory needs to be set, where the file "Sim1_Simulation.seed" and the script "Sim1_0_Functions_DataGeneration.R" for Simulation 1 or the file "Sim2_Simulation.seed" and the script "Sim2_0_Functions_DataGeneration.R" for Simulation 2 are saved. Furthermore, as mentioned the number of cores (num.cl) and number of iterations (n.sim) need to be kept in mind so that the running time is not increased too much. At last the folder path for saving of the generated data needs to be specified. 
  
- Data analysis of both simulation studies ("Sim1_2_DataAnalysis.R" and "Sim2_2_DataAnalysis.R"):\\
  These scripts conduct the data analysis of the generated data for all scenarios of both simulation studies. At the beginning the working directory needs to be set, where the script "Sim1_0_Functions_Analysis.R" for Simulation 1 or where the script "Sim2_0_Functions_Analysis.R" for Simulation 2 is saved. Furthermore, as mentioned the number of cores (num.cl) and number of iterations (n.sim) need to be kept in mind so that the running time is not increased too much. At last the folder path, where the generated data was saved (path_data), and the folder path, where the analysed data should be saved (path_ana) needs to be specified. 
  
- Optimal cutoff cutoff determination of both simulation studies ("Sim1_3_DataAnalysis_ASCO_CutoffValues.R" and "Sim2_3_DataAnalysis_ASCO_CutoffValues.R"): 
	These scripts conduct the investigating which ESMO/\IQWiG/\ModIQWiGspace category correspond to which ASCO score for both simulation studies. At the beginning the folder path, where the generated data and analysed data was saved, needs to be specified. As mentioned the number of cores (num.cl) and number of iterations (n.sim) need to be kept in mind so that the running time is not increased too much.
