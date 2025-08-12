###############################################################################
# SAS functions
###############################################################################
psas <- function(x, mu, sigma, epsilon, delta, log.p = FALSE){
  ifelse(sigma>0 & delta>0,
         logCDF <- pnorm(sinh(delta*asinh((x-mu)/sigma)-epsilon),0,1,log.p=TRUE),
         logCDF <- 'parameters out of range')
  ifelse( is.numeric(logCDF), ifelse(log.p, return(logCDF), return(exp(logCDF)) ), logCDF )
}

dsas <- function(x, mu, sigma, epsilon, delta, log=FALSE){
  ifelse(sigma>0 & delta>0,
         logPDF <- dnorm(sinh(delta*asinh((x-mu)/sigma)-epsilon),0,1,log=TRUE) +
           log(delta) + log(cosh(delta*asinh((x-mu)/sigma)-epsilon)) -
           0.5*log(1+(x-mu)^2/sigma^2) - log(sigma),
         logPDF <- 'parameters out of range')
  ifelse( is.numeric(logPDF), ifelse(log, return(logPDF), return(exp(logPDF)) ), logPDF )
}

qsas <- function(p, mu, sigma, epsilon, delta){
  ifelse(sigma>0 & delta>0,
         Q <- mu + sigma*sinh((asinh(qnorm(p)) + epsilon)/delta),
         Q <- 'parameters out of range')
  return(Q)
}

###############################################################################
# Monte Carlo ELBO
# phi layout: [1:4]=mu, [5:8]=log_sigma, [9:12]=epsilon, [13:16]=log_delta
###############################################################################
mc_elbo_SAS <- function(phi, N = 10000) {
  mu        <- phi[1:4]
  sigma     <- exp(phi[5:8])
  epsilon   <- phi[9:12]
  delta     <- exp(phi[13:16])
  
  U <- matrix(runif(N * 4), ncol = 4)
  
  eta <- cbind(
    qsas(U[,1], mu[1], sigma[1], epsilon[1], delta[1]),
    qsas(U[,2], mu[2], sigma[2], epsilon[2], delta[2]),
    qsas(U[,3], mu[3], sigma[3], epsilon[3], delta[3]),
    qsas(U[,4], mu[4], sigma[4], epsilon[4], delta[4])
  )
  
  # Term 1: E_q[log p(y, eta)]
  log_p <- numeric(N)
  for (i in 1:N) log_p[i] <- -log_postHRL(eta[i,])
  
  # Term 2: E_q[log q(eta)]
  log_q <- dsas(eta[,1], mu[1], sigma[1], epsilon[1], delta[1], log = TRUE) +
    dsas(eta[,2], mu[2], sigma[2], epsilon[2], delta[2], log = TRUE) +
    dsas(eta[,3], mu[3], sigma[3], epsilon[3], delta[3], log = TRUE) +
    dsas(eta[,4], mu[4], sigma[4], epsilon[4], delta[4], log = TRUE)
  
  elbo <- mean(log_p - log_q)
  return(-elbo)
}


###############################################################################
# MLE under SAS likelihood
###############################################################################

sim <- chainHR[ind, ]  # matrix with 4 columns (eta- or theta-scale; same scale used in dsas)

# Negative log-likelihood for SAS
nll_sas <- function(par, X) {
  mu      <- par[1:4]
  sigma   <- exp(par[5:8])   
  epsilon <- par[9:12]
  delta   <- exp(par[13:16])
  
  logq <- 
    dsas(X[,1], mu[1], sigma[1], epsilon[1], delta[1], log = TRUE) +
    dsas(X[,2], mu[2], sigma[2], epsilon[2], delta[2], log = TRUE) +
    dsas(X[,3], mu[3], sigma[3], epsilon[3], delta[3], log = TRUE) +
    dsas(X[,4], mu[4], sigma[4], epsilon[4], delta[4], log = TRUE)

  return(-sum(logq))
}

start_par <- rep(0, 16)

OPT_mle <- nlminb(start = start_par,
                  objective = function(p) nll_sas(p, sim),
                  control   = list(iter.max = 10000, trace = 1))
print(OPT_mle)

# Extract MLEs in the same parameterisation used by VB
mu_mle <- OPT_mle$par[1:4]
log_sigma_mle <- OPT_mle$par[5:8]   
epsilon_mle <- OPT_mle$par[9:12]
log_delta_mle <- OPT_mle$par[13:16]  


phi_init_mle <- c(mu_mle, log_sigma_mle, epsilon_mle, log_delta_mle)


###############################################################################
# Optimisation using MLE as a starting point
###############################################################################

out_vb <- nlminb(phi_init_mle, mc_elbo_SAS_MF, N = 5000,
                 control = list(iter.max = 10000, trace = 1))
print(out_vb)


phi_hat     <- out_vb$par
mu_hat      <- phi_hat[1:4]
sigma_hat   <- exp(phi_hat[5:8])
epsilon_hat <- phi_hat[9:12]
delta_hat   <- exp(phi_hat[13:16])

# Draw posterior samples from the SAS VB (on eta), then map to theta = exp(eta)
set.seed(123)
N_samp <- 10000
U_s    <- matrix(runif(N_samp*4), ncol = 4)

eta_s <- cbind(
  qsas(U_s[,1], mu_hat[1], sigma_hat[1], epsilon_hat[1], delta_hat[1]),
  qsas(U_s[,2], mu_hat[2], sigma_hat[2], epsilon_hat[2], delta_hat[2]),
  qsas(U_s[,3], mu_hat[3], sigma_hat[3], epsilon_hat[3], delta_hat[3]),
  qsas(U_s[,4], mu_hat[4], sigma_hat[4], epsilon_hat[4], delta_hat[4])
)

theta_sas <- as.data.frame(exp(eta_s))
colnames(theta_sas) <- c("lambda","kappa","alpha","beta")


###############################################################################
# Graphs for comparison
###############################################################################

# Plots for SAS MLE

d_sas_theta <- function(theta, i) {
  # i in {1,2,3,4} for lambda,kappa,alpha,beta (on eta-scale parameters)
  mu_i      <- mu_mle[i]
  sigma_i   <- exp(log_sigma_mle[i])
  epsilon_i <- epsilon_mle[i]
  delta_i   <- exp(log_delta_mle[i])
  eta <- log(theta)
  dsas(eta, mu_i, sigma_i, epsilon_i, delta_i, log = FALSE) * (1 / theta)
}

# Vectorized wrappers for curve()
f_theta_1 <- function(x) d_sas_theta(x, 1)  # lambda
f_theta_2 <- function(x) d_sas_theta(x, 2)  # kappa
f_theta_3 <- function(x) d_sas_theta(x, 3)  # alpha
f_theta_4 <- function(x) d_sas_theta(x, 4)  # beta


# Plots for all paramters and comparisons:

# lambda
plot(density(exp(post_napp[,1])),
     main = expression(lambda ~ ": Comparisons between Normal, MCMC, VB~SAS"),
     lwd = 2, col = "blue",
     xlab = expression(lambda), ylab = "Density", xlim = c(0, 3))
lines(density(lambdapHR), lwd = 2, col = "red")
lines(density(theta_sas$lambda), lwd = 2, col = "orange")
curve(f_theta_1(x), from = 1e-6, to = 3, add = TRUE, lwd = 2, col = "darkgreen", n = 1000)
legend("topright", c("Normal", "MCMC", "VB SAS", "SAS MLE"),
       col = c("blue", "red", "orange", "darkgreen"), lwd = 2)

# kappa
plot(density(exp(post_napp[,2])),
     main = expression(kappa ~ ": Comparisons between Normal, MCMC, VB~SAS"),
     lwd = 2, col = "blue",
     xlab = expression(kappa), ylab = "Density",
     ylim = c(0, 40))
lines(density(kappapHR), lwd = 2, col = "red")
lines(density(theta_sas$kappa), lwd = 2, col = "orange")
curve(f_theta_2(x), from = 1e-6, to = 2, add = TRUE, lwd = 2, col = "darkgreen", n = 1000)
legend("topright", c("Normal", "MCMC", "VB SAS", "SAS MLE"),
       col = c("blue", "red", "orange", "darkgreen"), lwd = 2)

# alpha
plot(density(exp(post_napp[,3])),
     main = expression(alpha ~ ": Comparisons between Normal, MCMC, VB~SAS"),
     lwd = 2, col = "blue",
     xlab = expression(alpha), ylab = "Density",
     xlim = c(0, 11))
lines(density(alphapHR), lwd = 2, col = "red")
lines(density(theta_sas$alpha), lwd = 2, col = "orange")
curve(f_theta_3(x), from = 1e-6, to = 11, add = TRUE, lwd = 2, col = "darkgreen", n = 1000)
legend("topright", c("Normal", "MCMC", "VB SAS", "SAS MLE"),
       col = c("blue", "red", "orange", "darkgreen"), lwd = 2)

# beta
plot(density(exp(post_napp[,4])),
     main = expression(beta ~ ": Comparisons between Normal, MCMC, VB~SAS"),
     lwd = 2, col = "blue",
     xlab = expression(beta), ylab = "Density",
     xlim = c(0, 10))
lines(density(betapHR), lwd = 2, col = "red")
lines(density(theta_sas$beta), lwd = 2, col = "orange")
curve(f_theta_4(x), from = 1e-6, to = 10, add = TRUE, lwd = 2, col = "darkgreen", n = 1000)
legend("topright", c("Normal", "MCMC", "VB SAS", "SAS MLE"),
       col = c("blue", "red", "orange", "darkgreen"), lwd = 2)
