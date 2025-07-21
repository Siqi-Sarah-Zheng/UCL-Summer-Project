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
# Define the parameter grids and the log posterior 
#################################################################################


#par1_grid <- seq(from = 1, to = 5,  length.out = 10)
#par2_grid <- seq(from = 1, to = 5,  length.out = 10)
#par3_grid <- seq(from = 1, to = 5,  length.out = 10)
#par4_grid <- seq(from = 1, to = 5,  length.out = 10)

#par1_grid <- seq(from = 0.242, to = 0.800, length.out = 10) 
#par2_grid <- seq(from = -2.727, to = -1.845, length.out = 10)
#par3_grid <- seq(from = 1.639, to = 2.239, length.out = 10)
#par4_grid <- seq(from = 1.273, to = 1.915, length.out = 10)

#par1_grid <- seq(from = MAP[1] - 3*sigmas_from_norm_approx[1], to = MAP[1] + 3*sigmas_from_norm_approx[1], length.out = 10)
#par2_grid <- seq(from = MAP[2] - 3*sigmas_from_norm_approx[2], to = MAP[2] + 3*sigmas_from_norm_approx[2], length.out = 10)
#par3_grid <- seq(from = MAP[3] - 3*sigmas_from_norm_approx[3], to = MAP[3] + 3*sigmas_from_norm_approx[3], length.out = 10)
#par4_grid <- seq(from = MAP[4] - 3*sigmas_from_norm_approx[4], to = MAP[4] + 3*sigmas_from_norm_approx[4], length.out = 10)

par1_grid <- seq(from = 0.2, to = 0.8,  length.out = 10)
par2_grid <- seq(from = -3, to = 3,  length.out = 10)
par3_grid <- seq(from = 1.5, to = 2.5,  length.out = 10)
par4_grid <- seq(from = 1, to = 2,  length.out = 10)

param_combinations <- expand.grid(par1 = par1_grid,
                                  par2 = par2_grid,
                                  par3 = par3_grid,
                                  par4 = par4_grid)
eta_grid  <- as.matrix(param_combinations)


N <- nrow(eta_grid)
log_p <- numeric(N)
for (i in 1:N) {
  log_p[i] <- -log_postHRL(eta_grid[i,])
}


#################################################################################
# LSE function
#################################################################################

LSE <- function(phi, eta, log_p){
  
  mu <- phi[1:4]
  sigma1 <- exp(phi[5:8])
  sigma2 <- exp(phi[9:12])
  
  if (!all(is.finite(sigma1)) || !all(is.finite(sigma2))) {
    return(.Machine$double.xmax)
  }
  
  # Calculate log_q using the pre-calculated grid.
  
  log_q_matrix <- matrix(NA, nrow = nrow(eta), ncol = 4)
  for (i in 1:4) {
    log_q_matrix[,i] <- dtp3(eta[,i], mu[i], sigma1[i], sigma2[i], FUN = dnorm, log = TRUE )
  }
  
  if (any(!is.finite(log_q_matrix))) {
    return(.Machine$double.xmax)
  }
  
  log_q <- rowSums(log_q_matrix)
  
  # Calculate MSE against the pre-calculated log_p vector.
  mse <- mean( (log_p - log_q)^2)
  
  return(mse)
}


#################################################################################
# Optimisation
#################################################################################

# phi_init <- c(as.numeric(inits), rep(0,8))

# We define the initial phi using the means and sd from the normal approximation, we set sigma 1 = sigma 2
sigmas_napp = sqrt(diag(Sigma))
phi_init_lse <- c(MAP, log(sigmas_napp), log(sigmas_from_napp))

set.seed(1234)

output_lse <- nlminb(phi_init, LSE, 
                 eta = eta_grid, 
                 log_p = log_p,
                 control = list(iter.max = 1e4, trace = 1))

phi_lse <- output$par

#################################################################################
# Result
#################################################################################

mu_lse = phi[1:4]
sigma1_lse = exp(phi[5:8])
sigma2_lse = exp(phi[9:12])

eta_samples_lse = cbind(rtp3(10000, mu_lse[1], sigma1_lse[1], sigma2_lse[1], FUN = rnorm),
                    rtp3(10000, mu_lse[2], sigma1_lse[2], sigma2_lse[2], FUN = rnorm),
                    rtp3(10000, mu_lse[3], sigma1_lse[3], sigma2_lse[3], FUN = rnorm),
                    rtp3(10000, mu_lse[4], sigma1_lse[4], sigma2_lse[4], FUN = rnorm))

theta_samples_lse <- as.data.frame(exp(eta_samples_lse))
colnames(theta_samples_lse) <- c("lambda", "kappa", "alpha", "beta")

summary(theta_samples_lse)



###################################################################################
# Graphs Comparison
###################################################################################


# Plot for lambda
plot(density(theta_samples_lse$lambda), main= expression(lambda ~ ": Normal, MCMC, VB , LSE"),
     lwd = 2, col = "darkgreen", xlab=expression(lambda), ylab = "Density", xlim = c(1,3))

lines(density(lambdapHR), lwd = 2, col = "red")
lines(density(theta_samples$lambda), lwd = 2, col = "purple")
lines(density(exp(post_napp[,1])), lwd = 2, col = "blue")
legend("topright", c("Normal","MCMC", "VB", "LSE"), col=c("blue","red", "purple", "darkgreen"), lwd=2)


# kappa
plot(density(exp(post_napp[,2])), main = expression(kappa ~ ": Normal, MCMC, VB , LSE"),
     lwd = 2, col = "blue", xlab = expression(kappa), ylab = "Density")
lines(density(kappapHR), lwd = 2, col = "red")
lines(density(theta_samples$kappa), lwd = 2, col = "purple")
lines(density(theta_samples_lse$kappa), lwd = 2, col = "darkgreen")
legend("topright", c("Normal","MCMC", "VB", "LSE"), col=c("blue","red", "purple", "darkgreen"), lwd=2)

# alpha
plot(density(exp(post_napp[,3])), main = expression(alpha ~ ": Normal, MCMC, VB , LSE"),
     lwd = 2, col = "blue", xlab = expression(alpha), ylab = "Density", xlim = c(0,11))
lines(density(alphapHR), lwd = 2, col = "red")
lines(density(theta_samples$alpha), lwd = 2, col = "purple")
lines(density(theta_samples_lse$alpha), lwd = 2, col = "darkgreen")
legend("topright", c("Normal","MCMC", "VB", "LSE"), col=c("blue","red", "purple", "darkgreen"), lwd=2)

# beta
plot(density(exp(post_napp[,4])), main = expression(beta ~ ": Normal, MCMC, VB , LSE"),
     lwd = 2, col = "blue", xlab = expression(beta), ylab = "Density", xlim = c(1,10), ylim = c(0,0.8))
lines(density(betapHR), lwd = 2, col = "red")
lines(density(theta_samples$beta), lwd = 2, col = "purple")
lines(density(theta_samples_lse$beta), lwd = 2, col = "darkgreen")
legend("topright", c("Normal","MCMC", "VB", "LSE"), col=c("blue","red", "purple", "darkgreen"), lwd=2)


