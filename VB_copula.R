library(mvtnorm)
# install.packages(c("deSolve", "survival", "devtools"))
library(deSolve)
library(survival)
library(devtools)
# install_github("FJRubio67/twopiece")
library(twopiece)
library(optimx)

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



# I adjusted the original qtp3 function to remove the NaNs in optimisation
qtp3_new <- function(p, mu, par1, par2, FUN, param = "tp") {

  clip_prob <- function(x, eps = 1e-10) {
    x[x < eps]     <- eps
    x[x > 1 - eps] <- 1 - eps
    x
  }

  param <- match.arg(param, choices = c("tp", "eps", "isf"))

  if (param == "tp") {
    if (par1 > 0 && par2 > 0) {
      Q <- ifelse(
        p < par1 / (par1 + par2),
        mu + par1 * FUN(clip_prob(0.5 * p * (par1 + par2) / par1)),
        mu + par2 * FUN(clip_prob(0.5 * ((par1 + par2) * (1 + p) - 2 * par1) / par2))
      )
    } else {
      stop("qtp3_new: par1 and par2 must both be > 0 for 'tp'")
    }
  }

  if (param == "eps") {
    sigma <- par1
    gamma <- par2
    if (sigma > 0 && abs(gamma) < 1) {
      Q <- ifelse(
        p < 0.5 * (1 + gamma),
        mu + sigma * (1 + gamma) * FUN(clip_prob(p / (1 + gamma))),
        mu + sigma * (1 - gamma) * FUN(clip_prob((p - gamma) / (1 - gamma)))
      )
    } else {
      stop("qtp3_new: sigma > 0 and |gamma| < 1 are required for 'eps'")
    }
  }

  if (param == "isf") {
    sigma <- par1
    gamma <- par2
    if (sigma > 0 && gamma > 0) {
      Q <- ifelse(
        p < gamma^2 / (1 + gamma^2),
        mu + sigma * gamma * FUN(clip_prob(0.5 * p * (1 + gamma^2) / gamma^2)),
        mu + sigma       * FUN(clip_prob(0.5 * (p * (1 + gamma^2) + 1 - gamma^2))) / gamma
      )
    } else {
      stop("qtp3_new: sigma > 0 and gamma > 0 are required for 'isf'")
    }
  }

  return(Q)
}



#################################################################################
# Monte Carlo function of the ELBO
#################################################################################


mc_elbo <- function(phi, N = 10000){
  
  # define the variational parameters
  mu = phi[1:4]
  sigma1 = exp(phi[5:8])
  sigma2 = exp(phi[9:12])
  psi <- phi[13:18]
  rhos <- 2 * plogis(psi) - 1
  
  R <- matrix(c(1, rhos[1], rhos[2], rhos[3],
                rhos[1], 1, rhos[4], rhos[5],
                rhos[2], rhos[4], 1, rhos[6],
                rhos[3], rhos[5], rhos[6], 1), nrow = 4, byrow = TRUE)

  Z1 <- mvrnorm(n = N/2, mu = rep(0, 4), Sigma = R)
  Z2 <- -Z1
  Z <- rbind(Z1, Z2)
  
  # Z <- mvrnorm(n = N, mu = rep(0, 4), Sigma = R)
  
  U <- pnorm(Z)
 
  eta <- cbind(
    qtp3_new(U[,1], mu[1], sigma1[1], sigma2[1], FUN = qnorm),
    qtp3_new(U[,2], mu[2], sigma1[2], sigma2[2], FUN = qnorm),
    qtp3_new(U[,3], mu[3], sigma1[3], sigma2[3], FUN = qnorm),
    qtp3_new(U[,4], mu[4], sigma1[4], sigma2[4], FUN = qnorm)
  )
  
  
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

phi_init = c(OPT$par, rep(0, 6))

set.seed(42)

output_coarse = nlminb(phi_init, mc_elbo, N = 1000, control = list(iter.max = 1e4, trace = 1)) # use N=100000

phi_init_fine = output_coarse$par

output = nlminb(phi_init_fine, mc_elbo, N=5000, control = list(iter.max = 1e4, trace = 1)) # use N=100000

# output = optimx(par = phi_init, fn = mc_elbo, N = 10000, method = "BFGS", control = list(trace = 1))

phi = output$par

print(phi)


#################################################################################
# Results
#################################################################################

#  parameters of the variational distribution
mu <- phi[1:4]
sigma1 <- exp(phi[5:8])
sigma2 <- exp(phi[9:12])
psi <- phi[13:18]
rhos <- 2 * plogis(psi) - 1

# Reconstruct the correlation matrix R
R <- matrix(c(1,      rhos[1], rhos[2], rhos[3],
              rhos[1], 1,       rhos[4], rhos[5],
              rhos[2], rhos[4], 1,       rhos[6],
              rhos[3], rhos[5], rhos[6], 1),
            nrow = 4, byrow = TRUE)


set.seed(42)

# Sample from the multivariate normal distribution of the copula
Z_samples <- mvrnorm(n = 10000, mu = rep(0, 4), Sigma = R)
U_samples <- pnorm(Z_samples)

eta_samples_copula <- cbind(
  qtp3_new(U_samples[,1], mu[1], sigma1[1], sigma2[1], FUN = qnorm),
  qtp3_new(U_samples[,2], mu[2], sigma1[2], sigma2[2], FUN = qnorm),
  qtp3_new(U_samples[,3], mu[3], sigma1[3], sigma2[3], FUN = qnorm),
  qtp3_new(U_samples[,4], mu[4], sigma1[4], sigma2[4], FUN = qnorm)
)


theta_copula_samples <- as.data.frame(exp(eta_samples_copula))
colnames(theta_copula_samples) <- c("lambda", "kappa", "alpha", "beta")

summary(theta_copula_samples)

#################################################################################
# Plots for Comparisons
#################################################################################

# Plot for lambda
plot(density(exp(post_napp[,1])),
     main = expression(lambda ~ ": Comparisons between Normal, MCMC, VB copula"),
     lwd = 2, col = "blue",
     xlab = expression(lambda), ylab = "Density")
lines(density(lambdapHR), lwd = 2, col = "red")
lines(density(theta_copula_samples$lambda), lwd = 2, col = "orange")
legend("topright", c("Normal", "MCMC", "VB copula"),
       col = c("blue", "red", "orange"), lwd = 2)

# kappa
plot(density(exp(post_napp[,2])),
     main = expression(kappa ~ ": Comparisons between Normal, MCMC, VB copula"),
     lwd = 2, col = "blue",
     xlab = expression(kappa), ylab = "Density",
     ylim = c(0, 35))
lines(density(kappapHR), lwd = 2, col = "red")
lines(density(theta_copula_samples$kappa), lwd = 2, col = "orange")
legend("topright", c("Normal", "MCMC", "VB copula"),
       col = c("blue", "red", "orange"), lwd = 2)

# alpha
plot(density(exp(post_napp[,3])),
     main = expression(alpha ~ ": Comparisons between Normal, MCMC, VB copula"),
     lwd = 2, col = "blue",
     xlab = expression(alpha), ylab = "Density",
     xlim = c(0, 11))
lines(density(alphapHR), lwd = 2, col = "red")
lines(density(theta_copula_samples$alpha), lwd = 2, col = "orange")
legend("topright", c("Normal", "MCMC", "VB copula"),
       col = c("blue", "red", "orange"), lwd = 2)

# beta
plot(density(exp(post_napp[,4])),
     main = expression(beta ~ ": Comparisons between Normal, MCMC, VB copula"),
     lwd = 2, col = "blue",
     xlab = expression(beta), ylab = "Density",
     xlim = c(1, 10))
lines(density(betapHR), lwd = 2, col = "red")
lines(density(theta_copula_samples$beta), lwd = 2, col = "orange")
legend("topright", c("Normal", "MCMC", "VB copula"),
       col = c("blue", "red", "orange"), lwd = 2)

