#!/usr/bin/env Rscript
# Author: Katie Bickerton kb620@kent.ac.uk
# Script: SIM_generate_ch.R
# Desc: Generate simulated data and run mark-recapture models with simulated data.
# Arguments: 
# Date: March 2022


# clear workspace
rm(list=ls())
graphics.off()


require(tidyverse)
require(lubridate)
require(gtools)

source("all_functions.R")


# 1. Generate simulated data ####

# JS translocation and JS standard run in Matlab
# Same capture histories used for all models


# Create simulated data


sim_data <- data.frame(matrix(ncol = 13, nrow = 0))


# scenario 1
sim_ch_S1 <- data.frame(matrix(ncol = 2, nrow = 0))
for(x in 1:250){
  
  caphist <- generate_capture_histories(N = 500, K = 10, tau = 4, p_int1 = 0.1, p_int2 = 0.3, phi_int1 = 0.476, phi_int2 = 0.886,
                                        N.tran = 15, t_int = rep(0.5,10), beta_int1 = 0.1, beta_int2 = 0.2)
  
  sim_totalch <- c("total_ch", caphist$total_ch, x, "S1")
  sim_totalpresent <- c("total_present", caphist$total_present, x, "S1")
  sim_betas <- c("beta", caphist$betas, x, "S1")
  sim_phi <- c("phi", caphist$phi, x, "S1")
  sim_p <- c("p", caphist$p, x, "S1")
  
  sim_data <- rbind(sim_data, as.data.frame(rbind(sim_totalch, sim_totalpresent, sim_betas, sim_phi, sim_p), row.names = FALSE))
  
  caphist$ch$sim <- rep(x, nrow(caphist$ch))
  caphist$ch <- caphist$ch[order(as.numeric(caphist$ch$ch), decreasing = TRUE),]
  sim_ch_S1 <- rbind(sim_ch_S1, caphist$ch)
  
  # write.csv(caphist$matrix_ch[1:15,], paste0("../Matlab/SIM/s1_simT", x, ".csv"), row.names = FALSE)
  # write.csv(caphist$matrix_ch[16:nrow(caphist$matrix_ch),], paste0("../Matlab/SIM/s1_simW", x, ".csv"), row.names = FALSE)
  
}

names(sim_ch_S1) <- c("ch","sim")


# scenario 2
sim_ch_S2 <- data.frame(matrix(ncol = 2, nrow = 0))

for(x in 1:250){
  
  caphist <- generate_capture_histories(N = 500, K = 10, tau = 4, p_int1 = 0.4, p_int2 = 0.6, phi_int1 = 0.476, phi_int2 = 0.886,
                                        N.tran = 15, t_int = rep(0.5,10), beta_int1 = 0.1, beta_int2 = 0.2)
  
  sim_totalch <- c("total_ch", caphist$total_ch, x, "S2")
  sim_totalpresent <- c("total_present", caphist$total_present, x, "S2")
  sim_betas <- c("beta", caphist$betas, x, "S2")
  sim_phi <- c("phi", caphist$phi, x, "S2")
  sim_p <- c("p", caphist$p, x, "S2")
  
  sim_data <- rbind(sim_data, as.data.frame(rbind(sim_totalch, sim_totalpresent, sim_betas, sim_phi, sim_p), row.names = FALSE))
  
  caphist$ch$sim <- rep(x, nrow(caphist$ch))
  caphist$ch <- caphist$ch[order(as.numeric(caphist$ch$ch), decreasing = TRUE),]
  sim_ch_S2 <- rbind(sim_ch_S2, caphist$ch)
  
  # write.csv(caphist$matrix_ch[1:15,], paste0("../Matlab/SIM/s2_simT", x, ".csv"), row.names = FALSE)
  # write.csv(caphist$matrix_ch[16:nrow(caphist$matrix_ch),], paste0("../Matlab/SIM/s2_simW", x, ".csv"), row.names = FALSE)
  
}

names(sim_ch_S2) <- c("ch","sim")

# scenario 3
sim_ch_S3 <- data.frame(matrix(ncol = 2, nrow = 0))

for(x in 1:250){
  
  caphist <- generate_capture_histories(N = 500, K = 10, tau = 4, p_int1 = 0.7, p_int2 = 0.9, phi_int1 = 0.476, phi_int2 = 0.886,
                                        N.tran = 15, t_int = rep(0.5,10), beta_int1 = 0.1, beta_int2 = 0.2)
  
  sim_totalch <- c("total_ch", caphist$total_ch, x, "S3")
  sim_totalpresent <- c("total_present", caphist$total_present, x, "S3")
  sim_betas <- c("beta", caphist$betas, x, "S3")
  sim_phi <- c("phi", caphist$phi, x, "S3")
  sim_p <- c("p", caphist$p, x, "S3")
  
  sim_data <- rbind(sim_data, as.data.frame(rbind(sim_totalch, sim_totalpresent, sim_betas, sim_phi, sim_p), row.names = FALSE))
  
  caphist$ch$sim <- rep(x, nrow(caphist$ch))
  caphist$ch <- caphist$ch[order(as.numeric(caphist$ch$ch), decreasing = TRUE),]
  sim_ch_S3 <- rbind(sim_ch_S3, caphist$ch)
  
  # write.csv(caphist$matrix_ch[1:15,], paste0("../Matlab/SIM/s3_simT", x, ".csv"), row.names = FALSE)
  # write.csv(caphist$matrix_ch[16:nrow(caphist$matrix_ch),], paste0("../Matlab/SIM/s3_simW", x, ".csv"), row.names = FALSE)
  
}

names(sim_ch_S3) <- c("ch","sim")


# scenario 4
sim_ch_S4 <- data.frame(matrix(ncol = 2, nrow = 0))

for(x in 1:250){
  
  caphist <- generate_capture_histories(N = 500, K = 10, tau = 4, p_int1 = 0.1, p_int2 = 0.3, phi_int1 = 0.476, phi_int2 = 0.886,
                                        N.tran = 30, t_int = rep(0.5,10), beta_int1 = 0.1, beta_int2 = 0.2)
  
  sim_totalch <- c("total_ch", caphist$total_ch, x, "S4")
  sim_totalpresent <- c("total_present", caphist$total_present, x, "S4")
  sim_betas <- c("beta", caphist$betas, x, "S4")
  sim_phi <- c("phi", caphist$phi, x, "S4")
  sim_p <- c("p", caphist$p, x, "S4")
  
  sim_data <- rbind(sim_data, as.data.frame(rbind(sim_totalch, sim_totalpresent, sim_betas, sim_phi, sim_p), row.names = FALSE))
  
  caphist$ch$sim <- rep(x, nrow(caphist$ch))
  caphist$ch <- caphist$ch[order(as.numeric(caphist$ch$ch), decreasing = TRUE),]
  sim_ch_S4 <- rbind(sim_ch_S4, caphist$ch)
  
  # write.csv(caphist$matrix_ch[1:30,], paste0("../Matlab/SIM/s4_simT", x, ".csv"), row.names = FALSE)
  # write.csv(caphist$matrix_ch[31:nrow(caphist$matrix_ch),], paste0("../Matlab/SIM/s4_simW", x, ".csv"), row.names = FALSE)
  
}

names(sim_ch_S4) <- c("ch","sim")


# scenario 5
sim_ch_S5 <- data.frame(matrix(ncol = 2, nrow = 0))

for(x in 1:250){
  
  caphist <- generate_capture_histories(N = 500, K = 10, tau = 4, p_int1 = 0.4, p_int2 = 0.6, phi_int1 = 0.476, phi_int2 = 0.886,
                                        N.tran = 30, t_int = rep(0.5,10), beta_int1 = 0.1, beta_int2 = 0.2)
  
  sim_totalch <- c("total_ch", caphist$total_ch, x, "S5")
  sim_totalpresent <- c("total_present", caphist$total_present, x, "S5")
  sim_betas <- c("beta", caphist$betas, x, "S5")
  sim_phi <- c("phi", caphist$phi, x, "S5")
  sim_p <- c("p", caphist$p, x, "S5")
  
  sim_data <- rbind(sim_data, as.data.frame(rbind(sim_totalch, sim_totalpresent, sim_betas, sim_phi, sim_p), row.names = FALSE))
  
  caphist$ch$sim <- rep(x, nrow(caphist$ch))
  caphist$ch <- caphist$ch[order(as.numeric(caphist$ch$ch), decreasing = TRUE),]
  sim_ch_S5 <- rbind(sim_ch_S5, caphist$ch)
  
  # write.csv(caphist$matrix_ch[1:30,], paste0("../Matlab/SIM/s5_simT", x, ".csv"), row.names = FALSE)
  # write.csv(caphist$matrix_ch[31:nrow(caphist$matrix_ch),], paste0("../Matlab/SIM/s5_simW", x, ".csv"), row.names = FALSE)
  
}

names(sim_ch_S5) <- c("ch","sim")


# scenario 6
sim_ch_S6 <- data.frame(matrix(ncol = 2, nrow = 0))

for(x in 1:250){
  
  caphist <- generate_capture_histories(N = 500, K = 10, tau = 4, p_int1 = 0.7, p_int2 = 0.9, phi_int1 = 0.476, phi_int2 = 0.886,
                                        N.tran = 30, t_int = rep(0.5,10), beta_int1 = 0.1, beta_int2 = 0.2)
  
  sim_totalch <- c("total_ch", caphist$total_ch, x, "S6")
  sim_totalpresent <- c("total_present", caphist$total_present, x, "S6")
  sim_betas <- c("beta", caphist$betas, x, "S6")
  sim_phi <- c("phi", caphist$phi, x, "S6")
  sim_p <- c("p", caphist$p, x, "S6")
  
  sim_data <- rbind(sim_data, as.data.frame(rbind(sim_totalch, sim_totalpresent, sim_betas, sim_phi, sim_p), row.names = FALSE))
  
  caphist$ch$sim <- rep(x, nrow(caphist$ch))
  caphist$ch <- caphist$ch[order(as.numeric(caphist$ch$ch), decreasing = TRUE),]
  sim_ch_S6 <- rbind(sim_ch_S6, caphist$ch)
  
  # write.csv(caphist$matrix_ch[1:30,], paste0("../Matlab/SIM/s6_simT", x, ".csv"), row.names = FALSE)
  # write.csv(caphist$matrix_ch[31:nrow(caphist$matrix_ch),], paste0("../Matlab/SIM/s6_simW", x, ".csv"), row.names = FALSE)
  
}

names(sim_ch_S6) <- c("ch","sim")


# scenario 7
sim_ch_S7 <- data.frame(matrix(ncol = 2, nrow = 0))

for(x in 1:250){
  
  caphist <- generate_capture_histories(N = 2000, K = 10, tau = 4, p_int1 = 0.1, p_int2 = 0.3, phi_int1 = 0.476, phi_int2 = 0.886,
                                        N.tran = 15, t_int = rep(0.5,10), beta_int1 = 0.1, beta_int2 = 0.2)
  
  sim_totalch <- c("total_ch", caphist$total_ch, x, "S7")
  sim_totalpresent <- c("total_present", caphist$total_present, x, "S7")
  sim_betas <- c("beta", caphist$betas, x, "S7")
  sim_phi <- c("phi", caphist$phi, x, "S7")
  sim_p <- c("p", caphist$p, x, "S7")
  
  sim_data <- rbind(sim_data, as.data.frame(rbind(sim_totalch, sim_totalpresent, sim_betas, sim_phi, sim_p), row.names = FALSE))
  
  caphist$ch$sim <- rep(x, nrow(caphist$ch))
  caphist$ch <- caphist$ch[order(as.numeric(caphist$ch$ch), decreasing = TRUE),]
  sim_ch_S7 <- rbind(sim_ch_S7, caphist$ch)
  
  # write.csv(caphist$matrix_ch[1:15,], paste0("../Matlab/SIM/s7_simT", x, ".csv"), row.names = FALSE)
  # write.csv(caphist$matrix_ch[16:nrow(caphist$matrix_ch),], paste0("../Matlab/SIM/s7_simW", x, ".csv"), row.names = FALSE)
  
}

names(sim_ch_S7) <- c("ch","sim")


# scenario 8
sim_ch_S8 <- data.frame(matrix(ncol = 2, nrow = 0))

for(x in 1:250){
  
  caphist <- generate_capture_histories(N = 2000, K = 10, tau = 4, p_int1 = 0.4, p_int2 = 0.6, phi_int1 = 0.476, phi_int2 = 0.886,
                                        N.tran = 15, t_int = rep(0.5,10), beta_int1 = 0.1, beta_int2 = 0.2)
  
  sim_totalch <- c("total_ch", caphist$total_ch, x, "S8")
  sim_totalpresent <- c("total_present", caphist$total_present, x, "S8")
  sim_betas <- c("beta", caphist$betas, x, "S8")
  sim_phi <- c("phi", caphist$phi, x, "S8")
  sim_p <- c("p", caphist$p, x, "S8")
  
  sim_data <- rbind(sim_data, as.data.frame(rbind(sim_totalch, sim_totalpresent, sim_betas, sim_phi, sim_p), row.names = FALSE))
  
  caphist$ch$sim <- rep(x, nrow(caphist$ch))
  caphist$ch <- caphist$ch[order(as.numeric(caphist$ch$ch), decreasing = TRUE),]
  sim_ch_S8 <- rbind(sim_ch_S8, caphist$ch)
  
  # write.csv(caphist$matrix_ch[1:15,], paste0("../Matlab/SIM/s8_simT", x, ".csv"), row.names = FALSE)
  # write.csv(caphist$matrix_ch[16:nrow(caphist$matrix_ch),], paste0("../Matlab/SIM/s8_simW", x, ".csv"), row.names = FALSE)
  
}

names(sim_ch_S8) <- c("ch","sim")


# scenario 9
sim_ch_S9 <- data.frame(matrix(ncol = 2, nrow = 0))

for(x in 1:250){
  
  caphist <- generate_capture_histories(N = 2000, K = 10, tau = 4, p_int1 = 0.7, p_int2 = 0.9, phi_int1 = 0.476, phi_int2 = 0.886,
                                        N.tran = 15, t_int = rep(0.5,10), beta_int1 = 0.1, beta_int2 = 0.2)
  
  sim_totalch <- c("total_ch", caphist$total_ch, x, "S9")
  sim_totalpresent <- c("total_present", caphist$total_present, x, "S9")
  sim_betas <- c("beta", caphist$betas, x, "S9")
  sim_phi <- c("phi", caphist$phi, x, "S9")
  sim_p <- c("p", caphist$p, x, "S9")
  
  sim_data <- rbind(sim_data, as.data.frame(rbind(sim_totalch, sim_totalpresent, sim_betas, sim_phi, sim_p), row.names = FALSE))
  
  caphist$ch$sim <- rep(x, nrow(caphist$ch))
  caphist$ch <- caphist$ch[order(as.numeric(caphist$ch$ch), decreasing = TRUE),]
  sim_ch_S9 <- rbind(sim_ch_S9, caphist$ch)
  
  # write.csv(caphist$matrix_ch[1:15,], paste0("../Matlab/SIM/s9_simT", x, ".csv"), row.names = FALSE)
  # write.csv(caphist$matrix_ch[16:nrow(caphist$matrix_ch),], paste0("../Matlab/SIM/s9_simW", x, ".csv"), row.names = FALSE)
  
}

names(sim_ch_S9) <- c("ch","sim")


# scenario 10
sim_ch_S10 <- data.frame(matrix(ncol = 2, nrow = 0))

for(x in 1:250){
  
  caphist <- generate_capture_histories(N = 2000, K = 10, tau = 4, p_int1 = 0.1, p_int2 = 0.3, phi_int1 = 0.476, phi_int2 = 0.886,
                                        N.tran = 30, t_int = rep(0.5,10), beta_int1 = 0.1, beta_int2 = 0.2)
  
  sim_totalch <- c("total_ch", caphist$total_ch, x, "S10")
  sim_totalpresent <- c("total_present", caphist$total_present, x, "S10")
  sim_betas <- c("beta", caphist$betas, x, "S10")
  sim_phi <- c("phi", caphist$phi, x, "S10")
  sim_p <- c("p", caphist$p, x, "S10")
  
  sim_data <- rbind(sim_data, as.data.frame(rbind(sim_totalch, sim_totalpresent, sim_betas, sim_phi, sim_p), row.names = FALSE))
  
  caphist$ch$sim <- rep(x, nrow(caphist$ch))
  caphist$ch <- caphist$ch[order(as.numeric(caphist$ch$ch), decreasing = TRUE),]
  sim_ch_S10 <- rbind(sim_ch_S10, caphist$ch)
  
  # write.csv(caphist$matrix_ch[1:30,], paste0("../Matlab/SIM/s10_simT", x, ".csv"), row.names = FALSE)
  # write.csv(caphist$matrix_ch[31:nrow(caphist$matrix_ch),], paste0("../Matlab/SIM/s10_simW", x, ".csv"), row.names = FALSE)
  
}

names(sim_ch_S10) <- c("ch","sim")


# scenario 11
sim_ch_S11 <- data.frame(matrix(ncol = 2, nrow = 0))

for(x in 1:250){
  
  caphist <- generate_capture_histories(N = 2000, K = 10, tau = 4, p_int1 = 0.4, p_int2 = 0.6, phi_int1 = 0.476, phi_int2 = 0.886,
                                        N.tran = 30, t_int = rep(0.5,10), beta_int1 = 0.1, beta_int2 = 0.2)
  
  sim_totalch <- c("total_ch", caphist$total_ch, x, "S11")
  sim_totalpresent <- c("total_present", caphist$total_present, x, "S11")
  sim_betas <- c("beta", caphist$betas, x, "S11")
  sim_phi <- c("phi", caphist$phi, x, "S11")
  sim_p <- c("p", caphist$p, x, "S11")
  
  sim_data <- rbind(sim_data, as.data.frame(rbind(sim_totalch, sim_totalpresent, sim_betas, sim_phi, sim_p), row.names = FALSE))
  
  caphist$ch$sim <- rep(x, nrow(caphist$ch))
  caphist$ch <- caphist$ch[order(as.numeric(caphist$ch$ch), decreasing = TRUE),]
  sim_ch_S11 <- rbind(sim_ch_S11, caphist$ch)
  
  # write.csv(caphist$matrix_ch[1:30,], paste0("../Matlab/SIM/s11_simT", x, ".csv"), row.names = FALSE)
  # write.csv(caphist$matrix_ch[31:nrow(caphist$matrix_ch),], paste0("../Matlab/SIM/s11_simW", x, ".csv"), row.names = FALSE)
  
}

names(sim_ch_S11) <- c("ch","sim")


# scenario 12
sim_ch_S12 <- data.frame(matrix(ncol = 2, nrow = 0))

for(x in 1:250){
  
  caphist <- generate_capture_histories(N = 2000, K = 10, tau = 4, p_int1 = 0.7, p_int2 = 0.9, phi_int1 = 0.476, phi_int2 = 0.886,
                                        N.tran = 30, t_int = rep(0.5,10), beta_int1 = 0.1, beta_int2 = 0.2)
  
  sim_totalch <- c("total_ch", caphist$total_ch, x, "S12")
  sim_totalpresent <- c("total_present", caphist$total_present, x, "S12")
  sim_betas <- c("beta", caphist$betas, x, "S12")
  sim_phi <- c("phi", caphist$phi, x, "S12")
  sim_p <- c("p", caphist$p, x, "S12")
  
  sim_data <- rbind(sim_data, as.data.frame(rbind(sim_totalch, sim_totalpresent, sim_betas, sim_phi, sim_p), row.names = FALSE))
  
  caphist$ch$sim <- rep(x, nrow(caphist$ch))
  caphist$ch <- caphist$ch[order(as.numeric(caphist$ch$ch), decreasing = TRUE),]
  sim_ch_S12 <- rbind(sim_ch_S12, caphist$ch)
  
  # write.csv(caphist$matrix_ch[1:30,], paste0("../Matlab/SIM/s12_simT", x, ".csv"), row.names = FALSE)
  # write.csv(caphist$matrix_ch[31:nrow(caphist$matrix_ch),], paste0("../Matlab/SIM/s12_simW", x, ".csv"), row.names = FALSE)
  
}

names(sim_ch_S12) <- c("ch","sim")
names(sim_data) <- c("Var",paste0("Time",seq(1, 10)), "sim","scenario")

save(sim_ch_S1, sim_ch_S2, sim_ch_S3, sim_ch_S4, sim_ch_S5, sim_ch_S6, sim_ch_S7, sim_ch_S8, sim_ch_S9, sim_ch_S10, 
     sim_ch_S11, sim_ch_S12, sim_data, file = "../Data/sim_ch.RData")




