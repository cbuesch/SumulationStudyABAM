# SumulationStudyABAM
Aus zweitem Paper --> noch anpassen!
Further information (R-Code, ADEMP structure of simulations) of Paper "A comparison of additional benefit assessment methods for time-to-event endpoints using hazard ratio point estimates or confidence interval limits by means of a simulation study".
In this github repository the paper and associated appendix including the ADEMP structure of our simulations can be found. Furthermore, the folder "FiguresAndTable" include all Figures and Tables of our paper as well as some plots illustrating further results of our simulations. And last but not least, the R-Code of our simulations are shared in the folder "RPrograms". Please read the R-Code instruction below before performing it!

## R-Code: 
Please keep in mind that the following R-packages should be installed: tidyverse, patchwork, ggpubr, survival, flexsurv, pcaPP, data.table, vcd, cutpointr, foreach and doParallel. The provided programs are using the doParallel package for parallel computing in windows systems. If you have a unix-like system you need to install doMC instead of doParallel and replace 
```r
cl <- makeCluster(num.cl)
registerDoParallel(cl)
```
by
```r
registerDoMC(num.cl)
```
Furthermore, you need to remove all lines containing 
```r
stopCluster(cl)
```
In addition, the running time of the programs (especially the data generation and data analysis) is very long. We performed it with parallel computing using 48 cores and still needed several days. Therefore, to reduce the running time, please use many cores for the parallel computing or reduce the number of iterations (n<sub>sim</sub>).

In the following several R-Code-Scripts are explained, which can be found in the folder "RPrograms":
- Costume functions ("0_CostumeFunctions_Analysis.R" and "0_CostumeFunctions_DataGeneration.R"):
  These two scripts contain all needed functions for the data generation and analysis. Hence, they need to be loaded before conducting the simulations. 
  
- Data Generation ("1_DataGeneration.R"):
  This script conducts the data generation for all four scenarios. At the beginning the working directory needs to be set, where the file "Simulation.seed" and the script "0_CostumeFunctions_DataGeneration.R" are saved. Furthermore, as mentioned above please keep the number of cores (num.cl) / number of iterations (n_sim) in mind so that the running time is not too long. At last the folder path for saving of the generated data should be supplied. 
  
- Data Analysis ("2_DataAnalysis.R"):
  This script conducts the data analysis of the generated data for all four scenarios. At the beginning the working directory needs to be set, where the script "0_CostumeFunctions_DataGeneration.R" is saved. Furthermore, as mentioned above please keep the number of cores (num.cl) in mind so that the running time is not too long. At last the folder path for saving of the analysed data (path_ana) and the folder path for the generated data (path_data) should be supplied. 
  
- Calculating ASCO cutoff values ("2_DataAnalysis_ASCO_CutoffValues.R"): This script conducts the investigating which ESMO / IQWiG<sub>RR</sub> / Mod-IQWiG<sub>HR</sub> category correspond to which ASCO score. At the beginning the working directory for the analysed data from "2_DataAnalysis.R" needs to be set. In addition, the folder path for saving of the calculated ASCO cutoff values (path_saving_cutoffs) need to be specified. As mentioned above please keep the number of cores (num.cl) in mind so that the running time is not too long.

- Visualizations ("3_VisualizationsOfResults.R"):
  This script creates all figures of the paper. In addition, it provides the Figures 1 and 2 of the Appendix. At the beginning you need to set three working directories. Firstly, to the folder where the final results of the script "2_DataAnalysis.R" is saved (path_ana). Secondly, to the folder where the generated data is saved (path_data). And thirdly, to the folder where the calculated ASCO cutoff values are saved (path_saving_cutoffs).
