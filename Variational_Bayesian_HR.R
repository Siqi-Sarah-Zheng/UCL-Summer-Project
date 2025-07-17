# install.packages(c("deSolve", "survival", "devtools"))
library(deSolve)
library(survival)
library(devtools)
# install_github("FJRubio67/twopiece")
library(twopiece)

source("routines.R")


#################################################################################
# Data preparation
#################################################################################

head(rotterdam)

dim(rotterdam)


# New data frame: logical status, time in years, survival times sorted
df <- data.frame(time = rotterdam$dtime, status = rotterdam$death)
df$status <- as.logical(rotterdam$death)
df$time <- df$time/365.24

df <- df[order(df$time),]

# Required quantities
status <- as.logical(df$status)
t_obs <- df$time[status]
survtimes <- df$time

#################################################################################
# Monte Carlo function of the ELBO
#################################################################################

mc_elbo <- function(phi, N = 1000){
  
  mu = phi[1:4]
  sigma1 = exp(phi[5:8])
  sigma2 = exp(phi[9:12])
  
  # Sample from the variational distribution
  
  # eps = matrix(rnorm(N*4), ncol = 4)
  
  eta = cbind(rtp3(N, mu[1], sigma1[1], sigma2[1], FUN = rnorm),
              rtp3(N, mu[2], sigma1[2], sigma2[2], FUN = rnorm),
              rtp3(N, mu[3], sigma1[3], sigma2[3], FUN = rnorm),
              rtp3(N, mu[4], sigma1[4], sigma2[4], FUN = rnorm))
  
  # Calculate the ELBO
  
  # Term 1: E[log p(y,θ)]
  
  log_p = numeric(N) # create an empty vector to store log p(y, θ)
  
  for (i in 1:N) {
    log_p[i] = -log_postHRL(eta[i,])
  }
  
  # Term 2: E[log q(η)]
  
  log_q_matrix = matrix(NA, nrow = N, ncol = 4)
  for (i in 1:4) {
    log_q_matrix[,i] = dtp3(eta[,i], mu[i], sigma1[i], sigma2[i], FUN = dnorm, log = TRUE )
  }
  log_q = rowSums(log_q_matrix)
  
  elbo = mean(log_p - log_q)
  
  return(-elbo)
}


#################################################################################
# Optimisation
#################################################################################
# phi_init = rep(0,12)
phi_init = c(as.numeric(inits), rep(0,8))

set.seed(42) # 42

phi = nlminb(phi_init, mc_elbo, control = list(iter.max = 1e4, trace = 1))$par

print(phi)


#################################################################################
# Results
#################################################################################

mu = phi[1:4]
sigma1 = exp(phi[5:8])
sigma2 = exp(phi[9:12])

eta_samples = cbind(rtp3(10000, mu[1], sigma1[1], sigma2[1], FUN = rnorm),
                    rtp3(10000, mu[2], sigma1[2], sigma2[2], FUN = rnorm),
                    rtp3(10000, mu[3], sigma1[3], sigma2[3], FUN = rnorm),
                    rtp3(10000, mu[4], sigma1[4], sigma2[4], FUN = rnorm))

theta_samples <- as.data.frame(exp(eta_samples))
colnames(theta_samples) <- c("lambda", "kappa", "alpha", "beta")

summary(theta_samples)
