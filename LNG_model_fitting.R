#!/usr/bin/env Rscript
# Author: Katie Bickerton kb620@kent.ac.uk
# Script: LNG_model_fitting.R
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

lng.ch.cov <- read.csv("IM_LNG_ch_covar.csv", header = TRUE, 
                       colClasses = c("ch" ="character"))
lng.ch.nocov <- read.csv("IM_LNG_ch_nocovar.csv", header = TRUE, 
                         colClasses = c("ch" ="character"))

lng.mat.cov <- read.csv("IM_LNG_mat_covar.csv", header = TRUE)
lng.mat.nocov <- read.csv("IM_LNG_mat_nocovar.csv", header = TRUE)

time_int <- read.csv("IM_LNG_time_intervals.csv", header = TRUE)
lng.env <- read.csv("IM_LNG_env_covars.csv", header = TRUE)


# Split capture histories and matrix capture histories by sex and whether they are translocated

lng.split <- sex.origin.split(lng.ch.cov)

lng.mat.split <- sex.origin.split(lng.mat.cov)

lng.ch.mf <- sex.split(lng.ch.cov)

# specify time intervals and time dependent covars
time.year <- head(time_int$year.int, -1)
airtemp <- lng.env$air.temp
moon <- lng.env$moon
effort <- time_int$effort.days
effort <- scale(effort,center=min(effort), scale=diff(range(effort)))[,1]

time.varying <- data.frame(effort, airtemp, moon)

## JS translocation model T
# input for translocation model is a list with standard chs (order MT,FT,MW,FW) 
# and matrix ch (order M,F)

# calculate number of parameters for model given specifications
n.theta <- calc.no.params(K = 20, tau = 5, time.covars = time.varying, time.beta = TRUE)
# theta with n values as 4 and remaining values to be estimated as 0
theta <- c(4,4, rep(0, (n.theta-2)))


y <- optim(par = theta, fn = js.translocation.run, ch.list = lng.split, mat.list = lng.mat.split, time.covars = time.varying,
           tau = 5, time_intervals = time.year, time.beta = TRUE, method = "BFGS", control = list(maxit = 10000))

## Standard JS model
# input split into male and female, standard ch and matrix formats

n.theta.JS <- calc.no.params(K = 20, tau = 1, time.covars = time.varying, time.beta = TRUE)
theta.JS <- c(4,4, rep(0, (n.theta.JS-2)))

w <- optim(par = theta.JS, fn = js.standard.run, ch.list = lng.ch.mf, mat.list = lng.mat.split, time.covars = time.varying,
           tau = 1, time_intervals = time.year, time.beta = TRUE, method = "BFGS", control = list(maxit = 10000))

## RMark models ####

require(RMark)

# run top model without covars
lng.rmark.1t <- select(lng.ch.cov, -ID,-sex,-origin)

lng.rmark.1tproc <- process.data(lng.rmark.1t, model = "POPAN", time.intervals = time.year)
lng.rmark.1tddl <- make.design.data(lng.rmark.1tproc)

rmark.1t <- mark(data = lng.rmark.1tproc, ddl = lng.rmark.1tddl, model = "POPAN", model.parameters = list(
  pent = list(formula=~time), p = list(formula=~1), Phi = list(formula=~time)))

rmark.1t.params <- rmark.param.calc1t(rmark.1t, K = 20)


# bootstrap confidence intervals

rmark.boot.params.1t <- list("beta" = data.frame(), "p" = data.frame(), "phi" = data.frame(),
                             "N" = data.frame(), "Nt" = data.frame())
names(rmark.boot.params.1t) <- c("beta","p","phi","N","Nt")


# run bootstrap for 250 runs
for(i in 1:250){
  
  boot.ch.1t <- generate.boot.ch(lng.rmark.1t)
  blng.rmark.1tproc <- process.data(boot.ch.1t, model = "POPAN", time.intervals = time.year)
  blng.rmark.1tddl <- make.design.data(blng.rmark.1tproc)
  
  rmark.1t.boot <- mark(data = blng.rmark.1tproc, ddl = blng.rmark.1tddl, model = "POPAN", model.parameters = list(
    pent = list(formula=~time), p = list(formula=~1), Phi = list(formula=~time)))
  
  brmark.1t.params <- rmark.param.calc1t(rmark.1t.boot, K = 20)
  
  rmark.boot.params.1t <- list.append.params(rmark.boot.params.1t, brmark.1t.params, i)
  
}

save(rmark.boot.params.1t,file =  "IM_lng_bootstrap_1t_rmark.RData")
save(rmark.1t.params, file = "IM_lng_parameters_1t_rmark.RData")


## marked models ####

## run model 1t

require(marked)

lng.marked <- lng.ch.cov %>% select(-ID, -sex, -origin)

marked.lng.proc <- process.data(lng.marked, model = "JS", time.intervals = time.year)
marked.lng.ddl <- make.design.data(marked.lng.proc)

marked.lng.1t <- crm(data = marked.lng.proc, ddl = marked.lng.ddl, model = "JS", time.intervals = time.year,
                     model.parameters = list(pent = list(formula=~time), p = list(formula=~1), 
                                             Phi = list(formula=~time)), hessian = FALSE,
                     accumulate = FALSE, itnmax = 1000)

marked.lng.params.1t <- calc.parameters.marked(marked.lng.1t)
names(marked.lng.params.1t) <- c("phi","p","beta","N","Nt")

marked.lng.params.1t$N$occ <- 1
marked.lng.params.1t$N$var.value <- "N"
marked.lng.params.1t$Nt <- data.frame("estimate" = marked.lng.params.1t$Nt, "occ" = seq(1,length(marked.lng.params.1t$Nt)))

marked.lng.params.1t$p$var.value <- "p"

marked.lng.params.1t$Nt$var.value <- rep("Nt",nrow(marked.lng.params.1t$Nt))
marked.lng.params.1t$beta$var.value <- rep("beta",nrow(marked.lng.params.1t$beta))
marked.lng.params.1t$phi$var.value <- rep("phi",nrow(marked.lng.params.1t$phi))

# run bootstrap

marked.boot.params.1t <- list("beta" = data.frame(), "p" = data.frame(), "Phi" = data.frame(),
                              "N" = data.frame(), "Nt" = data.frame())
names(marked.boot.params.1t) <- c("beta","p","Phi","N","Nt")


for(i in 1:250){
  
  boot.ch.1t.marked <- generate.boot.ch(lng.marked)
  blng.marked.1tproc <- process.data(boot.ch.1t.marked, model = "JS", time.intervals = time.year)
  blng.marked.1tddl <- make.design.data(blng.marked.1tproc)
  
  marked.1t.boot <- crm(data = blng.marked.1tproc, ddl = blng.marked.1tddl, model = "JS", time.intervals = time.year,
                        model.parameters = list(pent = list(formula=~time), p = list(formula=~1), 
                                                Phi = list(formula=~time)), hessian = FALSE,
                        accumulate = FALSE, itnmax = 1000)
  
  
  bmarked.1t.params <- calc.parameters.marked(marked.1t.boot)
  
  marked.boot.params.1t <- list.append.params.marked(marked.boot.params.1t, bmarked.1t.params, i)
  
}

names(marked.boot.params.1t)[[3]] <- "phi"

save(marked.boot.params.1t, file = "IM_lng_bootstrap_1t_marked.RData")
save(marked.lng.params.1t, file = "IM_lng_parameters_1t_marked.RData")







