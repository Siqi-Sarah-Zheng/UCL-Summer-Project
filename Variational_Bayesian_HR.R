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

mc_elbo <- function(phi, N = 5000){
  
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
  
  #log_p = numeric(N) # create an empty vector to store log p(y, θ)
  #for (i in 1:N) {
  #  log_p[i] = -log_postHRL(eta[i,])
  #}
  
  
  
  # Term 2: E[log q(η)]
  
  # log_q_matrix = matrix(NA, nrow = N, ncol = 4)
  
  log_q_matrix = matrix(NA, nrow = nrow(sim), ncol = 4)
  
  #for (i in 1:4) {
  #  log_q_matrix[,i] = dtp3(eta[,i], mu[i], sigma1[i], sigma2[i], FUN = dnorm, log = TRUE )
  #}
  
  for (j in 1:nrow(sim)) {
    for (i in 1:4) {
      log_q_matrix[j,i] <- dtp3(eta[j,i], mu[i], sigma1[i], sigma2[i], FUN = dnorm, log = TRUE )
    }
  }
  
  log_q = rowSums(log_q_matrix)
  
  elbo = mean(log_p - log_q)
  
  return(-elbo)
}


#################################################################################
# Optimisation
#################################################################################
# phi_init = rep(0,12)
# phi_init = c(as.numeric(inits), rep(0,8)) # mean and variance from normal approx. sigma1=sigma2=sigma from normal

sigmas_napp = sqrt(diag(Sigma))
phi_init <- c(MAP, log(sigmas_napp), log(sigmas_napp))


phi_init <- c(MAP, 
              log(c(sigmas_napp[1], 0.6 * sigmas_napp[2],
                    1.5 * sigmas_napp[3], 1.2 * sigmas_napp[4])), 
              log(c(1.5 * sigmas_napp[1], sigmas_napp[2], sigmas_napp[3], sigmas_napp[4])))


set.seed(42)

output = nlminb(phi_init, mc_elbo, control = list(iter.max = 1e4, trace = 1)) # use N=100000

phi = output$par

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


#################################################################################
# Plots for Comparisons
#################################################################################


# Plot for lambda
plot(density(exp(post_napp[,1])), main= expression(lambda ~ ": Comparisons between Normal, MCMC, VB"),
     lwd = 2, col = "blue", xlab=expression(lambda), ylab = "Density")
lines(density(lambdapHR), lwd = 2, col = "red")
lines(density(theta_samples$lambda), lwd = 2, col = "purple")
legend("topright", c("Normal","MCMC", "Variational Bayes"), col=c("blue","red", "purple"), lwd=2)


# kappa
plot(density(exp(post_napp[,2])), main = expression(kappa ~ ": Comparisons between Normal, MCMC, VB"),
     lwd = 2, col = "blue", xlab = expression(kappa), ylab = "Density", ylim = c(0,35))
lines(density(kappapHR), lwd = 2, col = "red")
lines(density(theta_samples$kappa), lwd = 2, col = "purple")
legend("topright", c("Normal","MCMC", "Variational Bayes"), col=c("blue","red", "purple"), lwd=2)


# alpha
plot(density(exp(post_napp[,3])), main = expression(alpha ~ ": Comparisons between Normal, MCMC, VB"),
     lwd = 2, col = "blue", xlab = expression(alpha), ylab = "Density", xlim = c(0,11))
lines(density(alphapHR), lwd = 2, col = "red")
lines(density(theta_samples$alpha), lwd = 2, col = "purple")
legend("topright", c("Normal","MCMC", "VB"), col=c("blue","red", "purple"), lwd=2) 


# beta
plot(density(exp(post_napp[,4])), main = expression(beta ~ ": Comparisons between Normal, MCMC, VB"),
     lwd = 2, col = "blue", xlab = expression(beta), ylab = "Density", xlim = c(1,10))
lines(density(betapHR), lwd = 2, col = "red")
lines(density(theta_samples$beta), lwd = 2, col = "purple")
legend("topright", c("Normal","MCMC", "VB"), col=c("blue","red", "purple"), lwd=2) 


