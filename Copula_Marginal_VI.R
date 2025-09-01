#################################################################################
# Libraries and Data preparation
#################################################################################
library(MASS)
library(mvtnorm)
library(deSolve)
library(survival)
library(twopiece)
library(numDeriv)
source("routines.R")

# rotterdam data
df <- data.frame(time = rotterdam$dtime, status = rotterdam$death)
df$status <- as.logical(rotterdam$death)
df$time   <- df$time/365.24
df <- df[order(df$time),]
status    <- as.logical(df$status)
t_obs     <- df$time[status]
survtimes <- df$time

load("Ex2.RData")   # Contains: inits, lambdapHR, kappapHR, alphapHR, betapHR

#################################################################################
# Normal approximation
#################################################################################
OPTMAP <- nlminb(inits, log_postHR, control = list(iter.max = 1e4))
MAP    <- OPTMAP$par
hess   <- hessian(func = log_postHR, x = MAP)
Sigma  <- solve(hess)

set.seed(123)
post_napp <- mvtnorm::rmvnorm(n = 1000, mean = MAP, sigma = Sigma)

#################################################################################
# Importance-Weighted Monte Carlo ELBO
#################################################################################
logmeanexp <- function(v){ m <- max(v); m + log(mean(exp(v - m))) }

mc_elbo <- function(phi, N = 5000, M = 10){
  N_use <- (N %/% M) * M
  mu   <- phi[1:4]
  s1   <- exp(phi[5:8])
  s2   <- exp(phi[9:12])
  psi  <- phi[13:18]
  rhos <- 2 * plogis(psi) - 1
  
  R <- matrix(c(
    1, rhos[1], rhos[2], rhos[3],
    rhos[1], 1, rhos[4], rhos[5],
    rhos[2], rhos[4], 1, rhos[6],
    rhos[3], rhos[5], rhos[6], 1
  ), 4, 4, byrow = TRUE)
  if (inherits(try(chol(R), silent = TRUE), "try-error")) return(.Machine$double.xmax)
  
  Z <- MASS::mvrnorm(n = N_use, mu = rep(0,4), Sigma = R)
  U <- pnorm(Z)
  
  eta <- cbind(
    twopiece::qtp3(U[,1], mu[1], s1[1], s2[1], FUN = qnorm),
    twopiece::qtp3(U[,2], mu[2], s1[2], s2[2], FUN = qnorm),
    twopiece::qtp3(U[,3], mu[3], s1[3], s2[3], FUN = qnorm),
    twopiece::qtp3(U[,4], mu[4], s1[4], s2[4], FUN = qnorm)
  )
  
  # E_q[log p(y,eta)]  (no offset)
  log_p <- vapply(seq_len(N_use), function(i) -log_postHRL(eta[i, ]), numeric(1))
  
  # E_q[log q(eta)] = ∑ log dtp3 + log copula
  log_q_marg <-
    twopiece::dtp3(eta[,1], mu[1], s1[1], s2[1], FUN = dnorm, log = TRUE) +
    twopiece::dtp3(eta[,2], mu[2], s1[2], s2[2], FUN = dnorm, log = TRUE) +
    twopiece::dtp3(eta[,3], mu[3], s1[3], s2[3], FUN = dnorm, log = TRUE) +
    twopiece::dtp3(eta[,4], mu[4], s1[4], s2[4], FUN = dnorm, log = TRUE)
  
  log_phi_R  <- mvtnorm::dmvnorm(Z, sigma = R, log = TRUE)
  log_phi_id <- rowSums(dnorm(Z, log = TRUE))
  log_q <- log_q_marg + (log_phi_R - log_phi_id)
  
  log_w <- log_p - log_q
  L <- matrix(log_w[1:((N_use/M)*M)], ncol = M, byrow = TRUE)
  iw_elbo <- mean(apply(L, 1, logmeanexp)) - log(M)
  -iw_elbo
}

#################################################################################
# Optimisation
#################################################################################
phi_init <- c(OPT$par, rep(0, 6))

set.seed(42)
out1 <- nlminb(phi_init, mc_elbo, N = 1000, M = 10, control = list(iter.max = 1e4, trace = 1))
out2 <- nlminb(out1$par, mc_elbo, N = 5000, M = 10, control = list(iter.max = 1e4, trace = 1))
phi  <- out2$par

#################################################################################
# Normalised Importance Sampling (NIS)
#################################################################################
mu     <- phi[1:4]
sigma1 <- exp(phi[5:8])
sigma2 <- exp(phi[9:12])
psi    <- phi[13:18]
rhos   <- 2 * plogis(psi) - 1

R <- matrix(c(1, rhos[1], rhos[2], rhos[3],
              rhos[1], 1, rhos[4], rhos[5],
              rhos[2], rhos[4], 1, rhos[6],
              rhos[3], rhos[5], rhos[6], 1), 4, 4, byrow = TRUE)

set.seed(123)
Zs <- mvrnorm(n = 10000, mu = rep(0, 4), Sigma = R)
Us <- pnorm(Zs)
etas <- cbind(
  twopiece::qtp3(Us[,1], mu[1], sigma1[1], sigma2[1], FUN = qnorm),
  twopiece::qtp3(Us[,2], mu[2], sigma1[2], sigma2[2], FUN = qnorm),
  twopiece::qtp3(Us[,3], mu[3], sigma1[3], sigma2[3], FUN = qnorm),
  twopiece::qtp3(Us[,4], mu[4], sigma1[4], sigma2[4], FUN = qnorm)
)

thetas <- as.data.frame(exp(etas))
names(thetas) <- c("lambda","kappa","alpha","beta")

# log q and log p
log_q_marg <- twopiece::dtp3(etas[,1], mu[1], sigma1[1], sigma2[1], FUN = dnorm, log = TRUE) +
  twopiece::dtp3(etas[,2], mu[2], sigma1[2], sigma2[2], FUN = dnorm, log = TRUE) +
  twopiece::dtp3(etas[,3], mu[3], sigma1[3], sigma2[3], FUN = dnorm, log = TRUE) +
  twopiece::dtp3(etas[,4], mu[4], sigma1[4], sigma2[4], FUN = dnorm, log = TRUE)
log_q_samp <- log_q_marg + (mvtnorm::dmvnorm(Zs, sigma = R, log = TRUE) - rowSums(dnorm(Zs, log = TRUE)))

log_p_samp <- vapply(seq_len(nrow(etas)), function(i) -log_postHRL(etas[i, ]), numeric(1))

# NIS weights
lw <- log_p_samp - log_q_samp
w  <- exp(lw - max(lw))
w  <- w / sum(w)

#################################################################################
# Graphs (Normal=blue, MCMC=red, VB=purple)
#################################################################################
par(mfrow=c(2,2), mar=c(4,4,2,1))

# λ
plot(density(exp(post_napp[,1])), main=expression(lambda), lwd=2, col="blue",
     xlab=expression(lambda), ylab="Density", xlim = range(c(lambdapHR, thetas$lambda)))
lines(density(lambdapHR), lwd=2, col="red")
lines(density(thetas$lambda, weights = w, bw = bw.nrd0(thetas$lambda)), lwd=2, col="purple")
legend("topright", c("Normal","MCMC","VB"), col=c("blue","red","purple"), lwd=2, bty="n", cex=0.85)

# κ
plot(density(exp(post_napp[,2])), main=expression(kappa), lwd=2, col="blue",
     xlab=expression(kappa), ylab="Density", ylim=c(0,26), xlim = range(c(kappapHR, thetas$kappa)))
lines(density(kappapHR), lwd=2, col="red")
lines(density(thetas$kappa, weights = w, bw = bw.nrd0(thetas$kappa)), lwd=2, col="purple")
legend("topright", c("Normal","MCMC","VB"), col=c("blue","red","purple"), lwd=2, bty="n", cex=0.85)

# α
plot(density(exp(post_napp[,3])), main=expression(alpha), lwd=2, col="blue",
     xlab=expression(alpha), ylab="Density", xlim = range(c(alphapHR, thetas$alpha)))
lines(density(alphapHR), lwd=2, col="red")
lines(density(thetas$alpha, weights = w, bw = bw.nrd0(thetas$alpha)), lwd=2, col="purple")
legend("topright", c("Normal","MCMC","VB"), col=c("blue","red","purple"), lwd=2, bty="n", cex=0.85)

# β
plot(density(exp(post_napp[,4])), main=expression(beta), lwd=2, col="blue",
     xlab=expression(beta), ylab="Density", xlim = range(c(betapHR, thetas$beta)))
lines(density(betapHR), lwd=2, col="red")
lines(density(thetas$beta, weights = w, bw = bw.nrd0(thetas$beta)), lwd = 2, col = "purple")
legend("topright", c("Normal","MCMC","VB"), col=c("blue","red","purple"), lwd=2, bty="n", cex=0.85)
