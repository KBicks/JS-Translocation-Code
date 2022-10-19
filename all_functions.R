#!/usr/bin/env Rscript
# Author: Katie Bickerton kb620@kent.ac.uk
# Script: all_functions.R
# Desc: Functions used in simulation study.
# Arguments: Null
# Date: May 2021




# generate simulated capture histories
generate_capture_histories <- function(N, K, tau, p_int1, p_int2, phi_int1, phi_int2, N.tran, t_int,
                                       beta_int1, beta_int2){
  
  # set empty matrices for capture history and presence histories
  PH <- matrix(0, N, K)
  CH <- matrix(0, N, K)
  
  # set entry probabilities for constant beta
  # betas <- c(rep(0, (tau-1)), 
  #             rep(1/(K-(tau-1)), K - (tau-1)))
  # 
  
  # set entry probabilities for time varying beta
  betas1 <- runif(K-(tau-1), min = beta_int1, max = beta_int2)
  betas2 <- betas1/sum(betas1)
  betas <- c(rep(0, tau-1), betas2)
  
  # set capture probabilities with stochasticity
  p <- runif(K-1, min = p_int1, max = p_int2)
  # set first capture to 1 for translocation release
  p <- c(1,p)
  
  # set phi probabilities with stochasticity
  phi <- runif(K, min = phi_int1 , max = phi_int2)
  phi <- phi^t_int
  
  
  
  # generate PH and CH for translocated individuals
  
  for(i in 1:N.tran){
    
    # entry time = 1
    en <- 1
    
    # generate exit time
    for(ex in en:K){
      if(rbinom(1, 1, phi[ex])==0) break
    }
    
    # generate presence history
    PH[i, en:ex] <- 1
    
    # generate capture history
    CH[i,] <- rbinom(K, 1, PH[i,]*p)
    
  }
  
  
  # generate PH and CH for non-translocated individuals
  for(i in (N.tran+1):N)
  {
    # generate entry time
    en <- which(rmultinom(1, 1, betas)==1)
    
    # generate exit time
    for(ex in en:K)
    {
      if(rbinom(1, 1, phi[ex])==0) break
    }
    
    # generate presence history (between entry and exit)
    PH[i, en:ex] <- 1 
    
    # generate capture history by multiplying presence by p
    CH[i,] <- rbinom(K, 1, PH[i,]*p)
    
  }
  
  # remove all 0 entries from CH
  CH <- CH[-which(apply(CH, 1, sum)==0),]
  matrix_CH <- CH
  
  # sample size
  # n <- nrow(CH)
  
  #check the number of translocations at time 1
  total_CH <- colSums(CH)
  
  
  # true number given perfect detection
  total_PH <- colSums(PH)
  
  CH <- data.frame(CH) %>% unite(col = "ch", 1:K, sep = "")
  outputs <- list(CH, total_CH, total_PH, betas, phi, p, matrix_CH)
  names(outputs) <- c("ch","total_ch","total_present", "betas","phi", "p", "matrix_ch")
  return(outputs)
  
}

# Run translocation JS model - required to be split by sex
js.translocation.run <- function(theta, ch.list, mat.list, time.covars, tau, time_intervals,
                                 sex.beta = FALSE, time.beta = FALSE, 
                                 time.p = FALSE, sex.p = FALSE, 
                                 time.phi = FALSE, sex.phi = FALSE){
  
  # load C++ likelihood using Rcpp package
  sourceCpp("JS_trans_llik.cpp")
  
  # Subsetting data for model ####
  
  mat_mt <- as.matrix(mat.list$MT)
  mat_ft <- as.matrix(mat.list$FT)
  mat_mw <- as.matrix(mat.list$MW)
  mat_fw <- as.matrix(mat.list$FW)
  
  
  K <- ncol(mat_mt)
  
  
  # If else statements so covariates take correct values ####
  
  # allow for sex dependent beta
  if(sex.beta == TRUE){
    sb <- 2
  }else{
    sb <- 1
  }
  
  # allow for time dependent beta
  if(time.beta == TRUE){
    tb <- K - tau
  }else{
    tb <- 1
  }
  
  # allow for sex p
  if(sex.p == TRUE){
    sp <- 2
  }else{
    sp <- 1
  }
  
  
  # allow for time dependent p
  if(time.p == TRUE){
    tp <- sp*(K-1)
    tpor <- 1
  }else{
    tp <- sp - 1
    tpor <- 0
  }
  
  # allow for time dependent covars for p
  if(length(time.covars > 0)){
    pcov <- sp*length(time.covars)
    mtp <- K-1
  }else{
    pcov <- sp - 1
    mtp <- 1
  }
  
  # allow for sex phi
  if(sex.phi == TRUE){
    sphi <- 2
  }else{
    sphi <- 1
  }
  
  
  # allow for time dependent phi
  if(time.phi == TRUE){
    tphi <- sphi*(K-1)
    mtphi <- K-1
  }else{
    tphi <- sphi - 1
    mtphi <- 1
  }
  
  
  
  
  # define theta positions for parameter optimisation ####
  theta_n <- theta[1:2] # n_m, n_f
  theta_b <- theta[3:((sb*tb)+2)] # betas
  theta_p <- theta[((sb*tb)+3):((sb*tb)+3+tp+pcov)] # p values with group and time covars
  theta_phi <- theta[((sb*tb)+4+tp+pcov):((sb*tb)+4+tp+pcov+tphi)] # phi values with group covars and time dependence
  
  
  # population size n ####
  nm <- exp(theta_n[1])       # n1 = N1-D1  never caught group 1 (what we want to estimate)
  Dm <- nrow(ch.list$MW)       #  Observed (wild) males
  Nm <- nm + Dm               #  N           total pop of caught+ uncaught in group 1
  
  nf <- exp(theta_n[2])       # n1 = N1-D1  never caught group 1 (what we want to estimate)
  Df <- nrow(ch.list$FW)       #  Observed (wild) females
  Nf <- nf + Df               #  N           total pop of caught+ uncaught in group 1
  
  
  # beta values ####
  theta_bM <- theta_b[1:(K-tau)]
  
  # accounting for whether sex effect present
  if(sex.beta == TRUE){
    theta_bF <- theta_b[(K-tau+1):length(theta_b)]
  }else{
    theta_bF <- theta_b[1:(K-tau)]
  }
  
  xpo_m <- exp(theta_bM) 
  xpo_f <- exp(theta_bF)
  
  beta_m <- c()
  beta_m[1:tau] <- 0                             
  beta_m[(tau+1):K] <- (xpo_m/(1+sum(xpo_m)))
  beta_m[tau] <- 1-sum(beta_m)
  
  beta_f <- c()
  beta_f[1:tau] <- 0
  beta_f[(tau+1):K] <- (xpo_f/(1+sum(xpo_f)))
  beta_f[tau] <- 1 - sum(beta_f)
  
  
  # p values ####
  
  # set initial values
  
  
  if(length(time.covars) > 0){
    
    p_mt <- 1
    p_ft <- 1
    p_mw <- rep(0, tau-1)
    p_fw <- rep(0, tau-1)
    
    if(length(time.covars)==1){
      
      
      if(sex.p == TRUE){
        
        for(i in 2:K){
          p_mt[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1])
          p_ft[i] <- inv.logit(theta_p[3]+theta_p[4]*time.covars[i,1])
        }
        
        p_mw[tau:K] <- p_mt[tau:K]
        p_fw[tau:K] <- p_ft[tau:K]
        
      }
      
      
      
      if(sex.p == FALSE){
        
        for(i in 2:K){
          p_mt[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1])
        }
        
        p_ft <- p_mt
        p_mw[tau:K] <- p_mt[tau:K]
        p_fw[tau:K] <- p_mt[tau:K]
        
      }
    }
    
    if(length(time.covars) == 2){
      
      
      if(sex.p == TRUE){
        
        for(i in 2:K){
          p_mt[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1]+theta_p[3]*time.covars[i,2])
          p_ft[i] <- inv.logit(theta_p[4]+theta_p[5]*time.covars[i,1]+theta_p[6]*time.covars[i,2])
        }
        p_mw[tau:K] <- p_mt[tau:K]
        p_fw[tau:K] <- p_ft[tau:K]
        
      }
      
      if(sex.p == FALSE){
        
        for(i in 2:K){
          p_mt[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1]+theta_p[3]*time.covars[i,2])
        }
        
        p_ft <- p_mt
        p_mw[tau:K] <- p_mt[tau:K]
        p_fw[tau:K] <- p_mt[tau:K]
        
      }
    }
    
    if(length(time.covars) == 3){
      
      if(sex.p == TRUE){
        
        for(i in 2:K){
          p_mt[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1]+theta_p[3]*time.covars[i,2]+theta_p[4]*time.covars[i,3])
          p_ft[i] <- inv.logit(theta_p[5]+theta_p[6]*time.covars[i,1]+theta_p[7]*time.covars[i,2]+theta_p[8]*time.covars[i,3])
        }
        p_mw[tau:K] <- p_mt[tau:K]
        p_fw[tau:K] <- p_ft[tau:K]
        
      }
      
      if(sex.p == FALSE){
        
        for(i in 2:K){
          p_mt[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1]+theta_p[3]*time.covars[i,2]+theta_p[4]*time.covars[i,3])
        }
        
        p_ft <- p_mt
        p_mw[tau:K] <- p_mt[tau:K]
        p_fw[tau:K] <- p_mt[tau:K]
        
      }
    }
    
    if(length(time.covars) == 4){
      
      if(sex.p == TRUE){
        
        for(i in 2:K){
          p_mt[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1]+theta_p[3]*time.covars[i,2]+theta_p[4]*time.covars[i,3]+theta_p[5]*time.covars[i,4])
          p_ft[i] <- inv.logit(theta_p[6]+theta_p[7]*time.covars[i,1]+theta_p[8]*time.covars[i,2]+theta_p[9]*time.covars[i,3]+theta_p[10]*time.covars[i,4])
        }
        p_mw[tau:K] <- p_mt[tau:K]
        p_fw[tau:K] <- p_ft[tau:K]
        
      }
      
      if(sex.p == FALSE){
        
        for(i in 2:K){
          p_mt[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1]+theta_p[3]*time.covars[i,2]+theta_p[4]*time.covars[i,3]+theta_p[5]*time.covars[i,4])
        }
        
        p_ft <- p_mt
        p_mw[tau:K] <- p_mt[tau:K]
        p_fw[tau:K] <- p_mt[tau:K]
        
      }
    }
  }else{
    p_mt <- c(1, inv.logit(theta_p[1:mtp]))
    
    if(sex.p == TRUE){
      p_ft <- c(1, inv.logit(theta_p[mtp+1:(2*mtp)]))
      p_mw <- p_mt
      p_fw <- p_ft
    }
    
    if(sex.p == FALSE){
      p_ft <- p_mt
      p_mw <- p_mt
      p_fw <- p_mt
    }
  }
  
  
  # phi values ####
  
  phi_mt <- c(inv.logit(theta_phi[1:mtphi])^time_intervals, 0)
  
  if(sex.phi == TRUE){
    phi_ft <- c(inv.logit(theta_phi[mtphi+1:(2*mtphi)])^time_intervals, 0)
    phi_mw <- phi_mt
    phi_fw <- phi_ft
  }
  
  if(sex.phi == FALSE){
    phi_ft <- phi_mt
    phi_mw <- phi_mt
    phi_fw <- phi_mt
  }
  
  
  # calculate entrance and exit points for males and females
  en_ex_mt <- calc.first.last(ch.list$MT)
  en_ex_ft <- calc.first.last(ch.list$FT)
  en_ex_mw <- calc.first.last(ch.list$MW)
  en_ex_fw <- calc.first.last(ch.list$FW)
  
  
  # Calculate likelihoods ####
  
  ### Likelihood of caught/uncaught wild males
  Lmt <- calc_llik_T(K=K, n=nm, N=Nm, en = en_ex_mt$en, ex = en_ex_mt$ex, 
                     p = p_mt, phi = phi_mt, matrix_ch = mat_mt)
  
  
  ### Likelihood of caught /uncaught wild females
  Lft <- calc_llik_T(K = K, n=nf, N=Nf, en = en_ex_ft$en, ex = en_ex_ft$ex, 
                     p = p_ft, phi = phi_ft, matrix_ch = mat_ft)
  
  
  ### Likelihood of caught/uncaught wild males
  Lmw <- calc_llik_W(tau=tau, K=K, n=nm, N=Nm, en = en_ex_mw$en, ex = en_ex_mw$ex, 
                     p = p_mw, phi = phi_mw, beta = beta_m, matrix_ch = mat_mw)
  
  
  ### Likelihood of caught /uncaught wild females
  Lfw <- calc_llik_W(tau=tau, K = K, n=nf, N=Nf, en = en_ex_fw$en, ex = en_ex_fw$ex, 
                     p = p_fw, phi = phi_fw, beta = beta_f, matrix_ch = mat_fw)
  
  
  ### Calculate full likelihood
  
  L <- Lmt + Lft + Lmw + Lfw
  return(L)
  
}

# Run standard JS model - required to be split by sex
js.standard.run <- function(theta, ch.list, mat.list, time.covars, tau, time_intervals,
                            sex.beta = FALSE, time.beta = FALSE, 
                            time.p = FALSE, sex.p = FALSE, 
                            time.phi = FALSE, sex.phi = FALSE){
  
  # load C++ likelihood using Rcpp package
  sourceCpp("JS_trans_llik.cpp")
  
  # Subsetting data for model ####
  
  mat_m <- as.matrix(rbind(mat.list$MT,mat.list$MW))
  mat_f <- as.matrix(rbind(mat.list$FT,mat.list$FW))
  
  
  K <- ncol(mat_m)
  
  
  # If else statements so covariates take correct values ####
  
  # allow for sex dependent beta
  if(sex.beta == TRUE){
    sb <- 2
  }else{
    sb <- 1
  }
  
  # allow for time dependent beta
  if(time.beta == TRUE){
    tb <- K - tau
  }else{
    tb <- 1
  }
  
  # allow for sex p
  if(sex.p == TRUE){
    sp <- 2
  }else{
    sp <- 1
  }
  
  
  # allow for time dependent p
  if(time.p == TRUE){
    tp <- sp*(K-1)
    tpor <- 1
  }else{
    tp <- sp - 1
    tpor <- 0
  }
  
  # allow for time dependent covars for p
  if(length(time.covars > 0)){
    pcov <- sp*length(time.covars)
    mtp <- K-1
  }else{
    pcov <- sp - 1
    mtp <- 1
  }
  
  # allow for sex phi
  if(sex.phi == TRUE){
    sphi <- 2
  }else{
    sphi <- 1
  }
  
  
  # allow for time dependent phi
  if(time.phi == TRUE){
    tphi <- sphi*(K-1)
    mtphi <- K-1
  }else{
    tphi <- sphi - 1
    mtphi <- 1
  }
  
  
  
  
  # define theta positions for parameter optimisation ####
  theta_n <- theta[1:2] # n_m, n_f
  theta_b <- theta[3:((sb*tb)+2)] # betas
  theta_p <- theta[((sb*tb)+3):((sb*tb)+3+tp+pcov)] # p values with group and time covars
  theta_phi <- theta[((sb*tb)+4+tp+pcov):((sb*tb)+4+tp+pcov+tphi)] # phi values with group covars and time dependence
  
  
  # population size n ####
  nm <- exp(theta_n[1])       # n1 = N1-D1  never caught group 1 (what we want to estimate)
  Dm <- nrow(ch.list$M)       #  Observed (wild) males
  Nm <- nm + Dm               #  N           total pop of caught+ uncaught in group 1
  
  nf <- exp(theta_n[2])       # n1 = N1-D1  never caught group 1 (what we want to estimate)
  Df <- nrow(ch.list$F)       #  Observed (wild) females
  Nf <- nf + Df               #  N           total pop of caught+ uncaught in group 1
  
  
  # beta values ####
  theta_bM <- theta_b[1:(K-tau)]
  
  # accounting for whether sex effect present
  if(sex.beta == TRUE){
    theta_bF <- theta_b[(K-tau+1):length(theta_b)]
  }else{
    theta_bF <- theta_b[1:(K-tau)]
  }
  
  xpo_m <- exp(theta_bM) 
  xpo_f <- exp(theta_bF)
  
  beta_m <- c()
  beta_m[1:tau] <- 0                             
  beta_m[(tau+1):K] <- (xpo_m/(1+sum(xpo_m)))
  beta_m[tau] <- 1-sum(beta_m)
  
  beta_f <- c()
  beta_f[1:tau] <- 0
  beta_f[(tau+1):K] <- (xpo_f/(1+sum(xpo_f)))
  beta_f[tau] <- 1 - sum(beta_f)
  
  
  # p values ####
  
  # set initial values
  
  
  if(length(time.covars) > 0){
    
    p_m <- c()
    p_f <- c()
    
    if(length(time.covars)==1){
      
      
      if(sex.p == TRUE){
        
        for(i in 1:K){
          p_m[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1])
          p_f[i] <- inv.logit(theta_p[3]+theta_p[4]*time.covars[i,1])
        }
        
      }
      
      
      
      if(sex.p == FALSE){
        
        for(i in 1:K){
          p_m[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1])
        }
        
        p_f <- p_m
        
      }
    }
    
    if(length(time.covars) == 2){
      
      
      if(sex.p == TRUE){
        
        for(i in 1:K){
          p_m[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1]+theta_p[3]*time.covars[i,2])
          p_f[i] <- inv.logit(theta_p[4]+theta_p[5]*time.covars[i,1]+theta_p[6]*time.covars[i,2])
        }
        
      }
      
      if(sex.p == FALSE){
        
        for(i in 1:K){
          p_m[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1]+theta_p[3]*time.covars[i,2])
        }
        
        p_f <- p_m
        
      }
    }
    
    if(length(time.covars) == 3){
      
      if(sex.p == TRUE){
        
        for(i in 1:K){
          p_m[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1]+theta_p[3]*time.covars[i,2]+theta_p[4]*time.covars[i,3])
          p_f[i] <- inv.logit(theta_p[5]+theta_p[6]*time.covars[i,1]+theta_p[7]*time.covars[i,2]+theta_p[8]*time.covars[i,3])
        }
        
      }
      
      if(sex.p == FALSE){
        
        for(i in 1:K){
          p_m[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1]+theta_p[3]*time.covars[i,2]+theta_p[4]*time.covars[i,3])
        }
        
        p_f <- p_m
        
      }
    }
    
    if(length(time.covars) == 4){
      
      if(sex.p == TRUE){
        
        for(i in 1:K){
          p_m[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1]+theta_p[3]*time.covars[i,2]+theta_p[4]*time.covars[i,3]+theta_p[5]*time.covars[i,4])
          p_f[i] <- inv.logit(theta_p[6]+theta_p[7]*time.covars[i,1]+theta_p[8]*time.covars[i,2]+theta_p[9]*time.covars[i,3]+theta_p[10]*time.covars[i,4])
        }
        
      }
      
      if(sex.p == FALSE){
        
        for(i in 1:K){
          p_m[i] <- inv.logit(theta_p[1]+theta_p[2]*time.covars[i,1]+theta_p[3]*time.covars[i,2]+theta_p[4]*time.covars[i,3]+theta_p[5]*time.covars[i,4])
        }
        
        p_f <- p_m
        
      }
    }
  }else{
    p_m <- c(1, inv.logit(theta_p[1:mtp]))
    
    if(sex.p == TRUE){
      p_f <- c(1, inv.logit(theta_p[mtp+1:(2*mtp)]))
    }
    
    if(sex.p == FALSE){
      p_f <- p_m
    }
  }
  
  
  # phi values ####
  
  phi_m <- c(inv.logit(theta_phi[1:mtphi])^time_intervals, 0)
  
  if(sex.phi == TRUE){
    phi_f <- c(inv.logit(theta_phi[mtphi+1:(2*mtphi)])^time_intervals, 0)
  }
  
  if(sex.phi == FALSE){
    phi_f <- phi_m
  }
  
  
  # calculate entrance and exit points for males and females
  en_ex_m <- calc.first.last(ch.list$M)
  en_ex_f <- calc.first.last(ch.list$F)
  
  
  # Calculate likelihoods ####
  
  ### Likelihood of caught/uncaught wild males
  Lm <- calc_llik_W(tau=tau, K=K, n=nm, N=Nm, en = en_ex_m$en, ex = en_ex_m$ex, 
                    p = p_m, phi = phi_m, beta = beta_m, matrix_ch = mat_m)
  
  
  ### Likelihood of caught /uncaught wild females
  Lf <- calc_llik_W(tau=tau, K = K, n=nf, N=Nf, en = en_ex_f$en, ex = en_ex_f$ex, 
                    p = p_f, phi = phi_f, beta = beta_f, matrix_ch = mat_f)
  
  
  ### Calculate full likelihood
  
  L <- Lm + Lf
  return(L)
  
  
  
  
}


# Run translocation JS model - not required to be split by sex, no covariates
# Run translocation JS model
js.translocation.run.nocov <- function(theta, ch.df, mat.df, tau, time_intervals,
                                 time.beta = FALSE,time.p = FALSE, 
                                 time.phi = FALSE){
  
  # load C++ likelihood using Rcpp package
  sourceCpp("JS_trans_llik.cpp")
  
  # Subsetting data for model ####
  
  mat_sim <- as.matrix(mat.df)
  
  K <- ncol(mat_sim)
  
  
  # If else statements for time dependence
  
  # allow for time dependent beta
  if(time.beta == TRUE){
    tb <- K - tau
  }else{
    tb <- 1
  }
  
  # allow for time dependent p
  if(time.p == TRUE){
    tp <- K-1
  }else{
    tp <- 1
  }
  
  
  # allow for time dependent phi
  if(time.phi == TRUE){
    tphi <- K-1
  }else{
    tphi <- 1
  }
  
  
  
  
  # define theta positions for parameter optimisation ####
  theta_n <- theta[1] # n
  theta_b <- theta[2:tb+2] # betas
  theta_p <- theta[tb+3:tb+tp+3] # p values with group and time covars
  theta_phi <- theta[tb+tp+4:tb+tp+tphi+4] # phi values with group covars and time dependence
  
  
  # population size n ####
  nm <- exp(theta_n[1])       # n1 = N1-D1  never caught group 1 (what we want to estimate)
  Dm <- nrow(ch.df)       #  Observed (wild) males
  Nm <- nm + Dm               #  N           total pop of caught+ uncaught in group 1
  
  
  # beta values ####
  theta_bM <- theta_b[1:(K-tau)]
  
  xpo_m <- exp(theta_bM) 
  
  beta_m <- c()
  beta_m[1:tau] <- 0                             
  beta_m[(tau+1):K] <- (xpo_m/(1+sum(xpo_m)))
  beta_m[tau] <- 1-sum(beta_m)
  
  # p values ####
  
  # set initial values
  
  p_m <- c(1, inv.logit(theta_p[1:tp]))
  
  # phi values ####
  
  phi_m <- c(inv.logit(theta_phi[1:tphi])^time_intervals, 0)
  
  
  # calculate entrance and exit points for males and females
  en_ex_m <- calc.first.last(ch.df)
  
  # Calculate likelihoods ####
  
  ### Likelihood of all
  L <- calc_llik_W(tau=tau, K = K, n=nm, N=Nm, en = en_ex_m$en, ex = en_ex_m$ex, 
                     p = p_m, phi = phi_m, beta = beta_m, matrix_ch = mat_sim)
  
  return(L)
  
}



sex.origin.split <- function(cmr.data){
  # split into separate data frames for sex and origin
  
  cmr.f <- split(cmr.data, cmr.data$sex)[[1]]
  cmr.m <- split(cmr.data, cmr.data$sex)[[2]]
  cmr.ft <- split(cmr.f, cmr.f$origin)[[1]] %>% select(-ID,-sex,-origin)
  cmr.fw <- split(cmr.f, cmr.f$origin)[[2]] %>% select(-ID,-sex,-origin)
  cmr.mt <- split(cmr.m, cmr.m$origin)[[1]] %>% select(-ID,-sex,-origin)
  cmr.mw <- split(cmr.m, cmr.m$origin)[[2]] %>% select(-ID,-sex,-origin)
  
  cmr.sep <- list(cmr.mt, cmr.ft, cmr.mw, cmr.fw)
  names(cmr.sep) <- c("MT","FT","MW","FW")
  return(cmr.sep)
  
}

sex.split <- function(cmr.data){
  # split into separate data frames for sex and origin
  
  cmr.f <- split(cmr.data, cmr.data$sex)[[1]] %>% select(-ID, -sex, -origin)
  cmr.m <- split(cmr.data, cmr.data$sex)[[2]] %>% select(-ID, -sex, -origin)
  
  cmr.sep <- list(cmr.m, cmr.f)
  names(cmr.sep) <- c("M","F")
  return(cmr.sep)
  
}


# calculate first and last sightings of each individual
calc.first.last <- function(ch.standard){
  
  en <- c()
  ex <- c()
  
  # for each individual
  for(i in 1:nrow(ch.standard)){
    
    # find entry and exit point
    en[i] <- min(unique(unlist(str_locate_all(ch.standard[i,], "1"))))
    ex[i] <- max(unique(unlist(str_locate_all(ch.standard[i,], "1"))))
  }
  
  en.ex <- list(en,ex)
  names(en.ex) <- c("en","ex")
  return(en.ex)
}

# calculate number of model parameters
calc.no.params <- function(K, tau, time.covars, sex.beta = FALSE, time.beta = FALSE, 
                           time.p = FALSE, sex.p = FALSE,
                           time.phi = FALSE, sex.phi = FALSE){
  
  # allow for sex dependent beta
  if(sex.beta == TRUE){
    sb <- 2
  }else{
    sb <- 1
  }
  
  # allow for time dependent beta
  if(time.beta == TRUE){
    tb <- K - tau
  }else{
    tb <- 1
  }
  
  # allow for sex p
  if(sex.p == TRUE){
    sp <- 2
  }else{
    sp <- 1
  }
  
  
  # allow for time dependent p
  if(time.p == TRUE){
    tp <- sp*(K-1)
  }else{
    tp <- sp - 1
  }
  
  # allow for time dependent covars for p
  if(length(time.covars > 0)){
    pcov <- sp*length(time.covars)
  }else{
    pcov <- sp - 1
  }
  
  # allow for sex phi
  if(sex.phi == TRUE){
    sphi <- 2
  }else{
    sphi <- 1
  }
  
  # allow for time dependent phi
  if(time.phi == TRUE){
    tphi <- sphi*(K-1)
    mtphi <- K-1
  }else{
    tphi <- sphi - 1
    mtphi <- 1
  }
  
  theta <- rep(0, ((sb*tb)+4+tp+pcov+tphi))
  return(length(theta))
  
}

rmark.param.calc1t <- function(model.output,K){
  
  # separate real param values
  params.rmark <- as.data.frame(model.output$results$real) %>% select(estimate,se,lcl,ucl)
  params.rmark$t <- c(seq(1,19),1,1,seq(2,20))
  params.rmark$var.value <- c(rep("phi",19),"p","N",rep("beta",19))
  beta.1 <- 1 - (sum(params.rmark$estimate[params.rmark$var.value == "beta"]))
  params.rmark[41,] <- c(beta.1, NA, NA, NA, 1, "beta")
  params.rmark$t <- as.numeric(params.rmark$t)
  params.rmark$estimate <- as.numeric(params.rmark$estimate)
  params.rmark <- params.rmark[order(params.rmark$t),]
  params.rmark <- list(params.rmark[params.rmark$var.value == "beta",],
                       params.rmark[params.rmark$var.value == "p",],
                       params.rmark[params.rmark$var.value == "phi",],
                       params.rmark[params.rmark$var.value == "N",])
  names(params.rmark) <- c("beta","p","phi","N")
  
  
  # calc Nt
  N1 <- params.rmark$N$estimate[1] * params.rmark$beta$estimate[1]
  Nt <- c(N1)
  
  for(i in 1:K-1){
    N <- (Nt[i]*(params.rmark$phi$estimate[i])) + (params.rmark$N$estimate[1] * params.rmark$beta$estimate[i+1])
    Nt <- append(Nt, N)
  }
  
  params.rmark$Nt <- data.frame("estimate" = Nt, "t" = seq(1,20), var.value = rep("Nt",K))
  
  return(params.rmark)
}


calc.parameters.marked <- function(models){
  
  # a <- models[[model.index]]
  a <- models
  model.params <- compute_real(a, se = TRUE)
  
  # calculate Nt
  pent0 <- 1 - sum(model.params$pent$estimate)
  
  N1 <- model.params$N$estimate[1] * pent0
  Nt <- c(N1)
  
  if(nrow(model.params$Phi) == 1){
    
    for(i in 1:19){
      N <- (Nt[i]*(model.params$Phi$estimate[1])) +  (model.params$N$estimate[1] * model.params$pent$estimate[i+1])
      Nt <- append(Nt, N)
    }
    
  }else{
    
    for(i in 1:19){
      N <- (Nt[i]*(model.params$Phi$estimate[i])) +  (model.params$N$estimate[1] * model.params$pent$estimate[i])
      Nt <- append(Nt, N)
    }
  }
  
  
  model.params$Nt <- Nt 
  
  return(model.params) 
}


generate.boot.ch <- function(data.ch){
  
  data.index <- floor(runif(nrow(data.ch), min = 1, max = nrow(data.ch))+0.5)
  boot.ch <- data.frame("ch" = as.character())
  
  for(i in 1:length(data.index)){
    
    boot.ch[i,1] <- data.ch[data.index[i],]
    
  }
  
  return(boot.ch)
  
}

list.append.params <- function(existing.params.list, new.params.list, index.no){
  
  new.params.list$beta <- new.params.list$beta %>% select(estimate, t, var.value)
  new.params.list$beta$run <- rep(index.no, nrow(new.params.list$beta))
  beta.frame <- rbind(existing.params.list$beta, new.params.list$beta)
  
  new.params.list$p <- new.params.list$p %>% select(estimate, t, var.value)
  new.params.list$p$run <- rep(index.no, nrow(new.params.list$p))
  p.frame <- rbind(existing.params.list$p, new.params.list$p)
  
  new.params.list$phi <- new.params.list$phi %>% select(estimate, t, var.value)
  new.params.list$phi$run <- rep(index.no, nrow(new.params.list$phi))
  phi.frame <- rbind(existing.params.list$phi, new.params.list$phi)
  
  new.params.list$N <- new.params.list$N %>% select(estimate, t, var.value)
  new.params.list$N$run <- rep(index.no, nrow(new.params.list$N))
  N.frame <- rbind(existing.params.list$N, new.params.list$N)
  
  new.params.list$Nt <- new.params.list$Nt %>% select(estimate, t, var.value)
  new.params.list$Nt$run <- rep(index.no, nrow(new.params.list$Nt))
  Nt.frame <- rbind(existing.params.list$Nt, new.params.list$Nt)
  
  existing.params.list <- list(beta.frame, p.frame, phi.frame, N.frame, Nt.frame)
  names(existing.params.list) <- c("beta","p","phi","N","Nt")
  
  return(existing.params.list)
}


list.append.params.marked <- function(existing.params.list, new.params.list, index.no){
  
  new.params.list$pent <- new.params.list$pent %>% select(estimate, time)
  new.params.list$pent$run <- rep(index.no, nrow(new.params.list$pent))
  new.params.list$pent$var.value <- rep("beta", nrow(new.params.list$pent))
  pent.frame <- rbind(existing.params.list$beta, new.params.list$pent)
  
  new.params.list$p <- new.params.list$p %>% select(estimate)
  new.params.list$p$run <- rep(index.no, nrow(new.params.list$p))
  new.params.list$p$var.value <- rep("p", nrow(new.params.list$p))
  p.frame <- rbind(existing.params.list$p, new.params.list$p)
  
  new.params.list$Phi <- new.params.list$Phi %>% select(estimate, time)
  new.params.list$Phi$run <- rep(index.no, nrow(new.params.list$Phi))
  new.params.list$Phi$var.value <- rep("phi", nrow(new.params.list$Phi))
  phi.frame <- rbind(existing.params.list$Phi, new.params.list$Phi)
  
  new.params.list$N <- new.params.list$N %>% select(estimate)
  new.params.list$N$run <- rep(index.no, nrow(new.params.list$N))
  new.params.list$N$var.value <- rep("N", nrow(new.params.list$N))
  N.frame <- rbind(existing.params.list$N, new.params.list$N)
  
  new.params.list$Nt <- data.frame("estimate" = new.params.list$Nt)
  new.params.list$time <- seq(1,20)
  new.params.list$Nt$run <- rep(index.no, nrow(new.params.list$Nt))
  new.params.list$Nt$var.value <- rep("Nt", nrow(new.params.list$Nt))
  Nt.frame <- rbind(existing.params.list$Nt, new.params.list$Nt)
  
  existing.params.list <- list(pent.frame, p.frame, phi.frame, N.frame, Nt.frame)
  names(existing.params.list) <- c("beta","p","Phi","N","Nt")
  
  return(existing.params.list)
}



