#!/usr/bin/env Rscript
# Author: Katie Bickerton kb620@kent.ac.uk
# Script: SIM_model_fitting.R
# Desc: Runs JS models and model output analysis of IM LNG data.
# Arguments: Capture histories file and covariates.
# Date: May 2021

# clear workspace
rm(list=ls())
graphics.off()


require(tidyverse)
require(lubridate)
require(gtools)
require(Rcpp)

source("all_functions.R")

load("sim_ch.RData")
sim_data[,2:11] <- lapply(sim_data[,2:11], as.numeric)


# see MATLAB code for model fitting (to be updated)


