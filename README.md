## Code files for *"Estimating of translocated populations: a modification of the Jolly-Seber model"*

**Authors:** Katherine T. Bickerton, John G. Ewen, Stefano Canessa, Nik C. Cole, Fay Frost, Rouben Mootoocurpen and Rachel McCrea

**Journal:** Ecological Applications

*Code File Descriptions:*  
**all_functions.R** - functions required for the running of all R scripts within this repository.  
**SIM_generate_ch.R** - generates capture histories for the 12 simulated scenarios in the simulation study and saves to an RData file. Options to save in format for MATLAB are available but commented out.  
**SIM_model_fitting.R** - fits standard POPAN formulation of Jolly-Seber model and modified translocation formulation of Jolly-Seber model to simulated data. This may take time to run due to optimisation in R.    
**LNG_model_fitting.R** - fits standard POPAN formulation of Jolly-Seber model and modified translocation formulation of Jolly-Seber model to case study data. Also fits models in packages *RMark* and *marked* with bootstrap confidence intervals. This may take some time to run due to optimisation in R.  
**JS_trans_llik.cpp** - C++ likelihood underlying the standard and translocation formulations of the Jolly-Seber model.  
**MATLAB** - Folder contains files to run models using MATLAB which is faster than using R. All files in the folder are either sourced from the **runLNG.m** file for the case study data or the **runSIMJS.m** file for the simulation study.  
**Data** - Folder contains csv files to run models using IM_LNG_run.R. All files in the folder are either sourced from the **IM_LNG_run.R** file for the case study, only processed data files are provided therefore data cleaning steps in **IM_LNG_run.R** are not required and users should begin on step 10.  


*Required R packages:*  
**tidyverse**  
**lubridate**  
**gtools**  
**Rcpp**  
**RMark**  
**marked**  

  [![DOI](https://zenodo.org/badge/553999125.svg)](https://zenodo.org/badge/latestdoi/553999125)

*Program Versions:*  
**R** version 4.1.2  
**MATLAB** version 9.8.0 (R2020a)



 
