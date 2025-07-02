library(torch)
library(deSolve)
library(survival)

source("routines.R")

###############################################################################
# Data preparation                                                            
###############################################################################

# New data frame: logical status, time in years, survival times sorted
df <- data.frame(time = rotterdam$dtime, status = rotterdam$death)
df$status <- as.logical(rotterdam$death)
df$time <- df$time/365.24
df <- df[order(df$time),]

# Required quantities
t_obs <- df$time                   
delta <- as.numeric(df$status)    
n     <- length(t_obs)
# 
# ###############################################################################
# # Hazard–Response ODE                                           
# ###############################################################################
# 
# # Hazard-Response ODE, from 'routines.R'
# hazmodHRL <- function(t, y, par) {
#   # state variables
#   lh <- y[1]
#   lq <- y[2]
#   CH <- y[3]
#   
#   # parameters
#   lambda <- par[1]
#   kappa <- par[2]
#   alpha <- par[3]
#   beta <- par[4]
#   
#   # model equations
#   dlh <-  lambda*(1 - exp(lh)/kappa) - alpha*exp(lq)
#   dlq <-  beta*(1-exp(lq)/kappa) - alpha*exp(lh)
#   dCH <- exp(lh)
#   
#   # result
#   return( list(c(dlh, dlq, dCH)) )
# }
# 
# 
# # Jacobian of the system of ODEs, from 'routines.R'
# 
# jacODE <- function(t, y, par) {
#   # state variables
#   h <- y[1]
#   q <- y[2]
#   CH <- y[3]
#   
#   # parameters
#   lambda <- par[1]
#   kappa <- par[2]
#   alpha <- par[3]
#   beta <- par[4]
#   
#   J1 <- c(lambda*(1-2*h/kappa) - alpha*q, -alpha*h, 0)
#   J2 <- c(-alpha*q, beta*(1-2*q/kappa) - alpha*h, 0)
#   J3 <- c(1, 0, 0)
#   
#   jacob <- as.matrix(rbind(J1, J2, J3))
#   return(jacob)
# }
# 
# jacODEL <- function(t, y, par) {
#   # state variables
#   lh <- y[1]
#   lq <- y[2]
#   CH <- y[3]
#   
#   # parameters
#   lambda <- par[1]
#   kappa <- par[2]
#   alpha <- par[3]
#   beta <- par[4]
#   
#   J1 <- c( - lambda*exp(lh)/kappa, - alpha*exp(lq), 0)
#   J2 <- c(- alpha*exp(lh),- beta*exp(lq)/kappa, 0)
#   J3 <- c(exp(lh), 0, 0)
#   
#   jacob <- as.matrix(rbind(J1, J2, J3))
#   return(jacob)
# }
# 
# 
# # Solve the Hazard-Response ODE system
# solve_HR <- function(theta, times, h0 = 0.1, q0 = 0.1) {
#   # Initial conditions in log scale
#   init <- c(
#     lh = log(h0),  # Initial log hazard
#     lq = log(q0),  # Initial log treatment response
#     CH = 0  # Cumulative hazard starts at 0
#   )
#   
#   # Solve ODE using LSODE method
#   out  <- ode(
#     y       = init,  # Initial state
#     times   = sort(unique(times)),  # Time points (include 0 for stability)
#     func    = hazmodHRL,
#     parms   = theta,  
#     jacfunc = jacODEL,
#     jactype = "fullusr",  
#     method  = "lsode",
#     rtol    = 1e-6,
#     atol    = 1e-6,
#     maxsteps= 5e4  # Maximum solver steps
#   )
#   
#   # Remove initial time point (0) if present
#   if (out[1,"time"] == 0) out <- out[-1, ]
#   
#   # Extract hazard and cumulative hazard
#   h <- pmax(exp(out[,"lh"]), 1e-8) 
#   H <- out[,"CH"]
#   
#   list(
#     h = approx(out[,"time"], h, xout = times)$y,  # Hazard at observation times
#     H = approx(out[,"time"], H, xout = times)$y  # Cumulative hazard at observation times
#   )
# }
# 
# 
# # Log-likelihood function
# loglik_theta <- function(theta) {
#   # Solve ODE system for current parameters
#   sol <- solve_HR(theta, t_obs)
#   
#   if (any(!is.finite(sol$h)) || any(!is.finite(sol$H))) {
#     return(-1e6)  # Return large negative value for invalid solutions
#   }
#   
#   # Calculate log-likelihood:
#   sum(delta * log(sol$h) - sol$H)
# }

###############################################################################
# Monte-Carlo ELBO                                   
###############################################################################

# ELBO calculation using torch
elbo_torch <- function(phi, S = 16) {
  
  # Variational parameters
  mu        <- phi$mu   # Log mean of variational distribution
  log_sd    <- phi$log_sd  # Log standard deviation of variational distribution
  sd        <- torch_exp(log_sd)  # Standard deviation
  
  # Reparameterisation trick: sample from variational distribution
  eps <- torch_randn(c(S, 4))  # Standard normal samples
  z <- mu$unsqueeze(1) + sd$unsqueeze(1) * eps  # Transform to log-parameter space
  theta <- torch_exp(z)  # Convert to original parameter space
  
  logp_vec <- torch_empty(S)  # log p(y, theta)
  logq_vec <- torch_empty(S)  # log q(theta; phi)
  
  # Calculate for each sample
  for (s in 1:S) {
    # Extract parameters for this sample
    theta_s <- as.numeric(theta[s, ])
    
    # Calculate log joint probability
    lp <- loglik_theta(theta_s) +
      sum(dgamma(theta_s, shape = 2, rate = 0.5, log = TRUE))
    logp_vec[s] <- lp
    
    # Calculate variational density: log q(theta; phi), using reparameterisation: z = mu + sd * eps
    logq_vec[s] <- torch_sum(
      -0.5 * ((z[s, ] - mu)/sd)**2 -
        log_sd -
        0.5 * log(2*pi)
    )
  }
  
  # ELBO = E_q[log p(y,theta) - log q(theta; phi)]
  torch_mean(logp_vec - logq_vec)
}

###############################################################################
# Optimisation                                               
###############################################################################
# Initialised variational parameters
phi <- list(
  mu     = torch_zeros(4, requires_grad = TRUE),
  log_sd = torch_full(c(4), -1, requires_grad = TRUE)
)

# Adam optimisor
opt <- optim_adam(list(phi$mu, phi$log_sd), lr = 2e-2)

# Optimisation loop
for (iter in 1:1000) {
  # Reset gradients
  opt$zero_grad()
  
  # Compute ELBO
  elbo <- elbo_torch(phi, S = 20)
  
  # Maximize ELBO is equivalent to minimize negative ELBO
  (-elbo)$backward()
  
  # Update parameters
  opt$step()
  
  # Print progress
  if (iter %% 100 == 0)
    cat(sprintf("iter %4d  ELBO %.3f\n", iter, elbo$item()))
}
