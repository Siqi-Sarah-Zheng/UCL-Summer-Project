rm(list=ls())
library(copula)
library(twopiece)
library(mvtnorm)


load("Ex2.RData")
sim <- chainHR[ind,]
#sim <- rmvnorm(n = 1000, mean = rep(2,4), sigma = diag(4))

# log likelihood function (reparameterised)
loglik = function(par){
  mu = par[1:4]
  sigma1 = exp(par[5:8])
  sigma2 = exp(par[9:12])
  psi <- par[13:18]
  rhos <- 2 * plogis(psi) - 1
  # Ingredients of the log likelihood
  norm.cop = normalCopula(rhos, dim = 4, dispstr = "un")
  probs =cbind(ptp3(sim[,1],mu[1],sigma1[1],sigma2[1],param="tp",FUN=pnorm),
               ptp3(sim[,2],mu[2],sigma1[2],sigma2[2],param="tp",FUN=pnorm),
               ptp3(sim[,3],mu[3],sigma1[3],sigma2[3],param="tp",FUN=pnorm),
               ptp3(sim[,4],mu[4],sigma1[4],sigma2[4],param="tp",FUN=pnorm))
  val.cop = sum(dCopula(probs, copula = norm.cop, log = TRUE))
  val.marg1 = sum( dtp3(sim[,1],mu[1],sigma1[1],sigma2[1],param="tp",FUN=dnorm, log = TRUE) )
  val.marg2 = sum( dtp3(sim[,2],mu[2],sigma1[2],sigma2[2],param="tp",FUN=dnorm, log = TRUE) )
  val.marg3 = sum( dtp3(sim[,3],mu[3],sigma1[3],sigma2[3],param="tp",FUN=dnorm, log = TRUE) )
  val.marg4 = sum( dtp3(sim[,4],mu[4],sigma1[4],sigma2[4],param="tp",FUN=dnorm, log = TRUE) )

  
  # output
  out <- -val.cop - val.marg1 - val.marg2 - val.marg3 - val.marg4
  
  return(out)
}


start = rep(0,18)
loglik(start)
# Optimisation step
OPT = nlminb(start,loglik,control=list(iter.max=1000))
OPT

MLE  <- c(OPT$par[1:4],exp(OPT$par[5:12]), 2*plogis(OPT$par[13:18])-1)

cbind(MLE, c(rep(2,4),rep(1,8),rep(0,6)))



# Simulate from the MLE estimates
mu     <- OPT$par[1:4]
sigma1 <- exp(OPT$par[5:8])
sigma2 <- exp(OPT$par[9:12])
psi    <- OPT$par[13:18]
rhos   <- 2*plogis(psi) - 1

R <- matrix(c(1,        rhos[1], rhos[2], rhos[3],
              rhos[1],  1,       rhos[4], rhos[5],
              rhos[2],  rhos[4], 1,       rhos[6],
              rhos[3],  rhos[5], rhos[6], 1),
            nrow = 4, byrow = TRUE)

Z_sample <- rmvnorm(5000, sigma = R)
U_sample <- pnorm(Z_sample)

eta_sample <- cbind(
  qtp3(U_sample[,1], mu[1], sigma1[1], sigma2[1], FUN = qnorm),
  qtp3(U_sample[,2], mu[2], sigma1[2], sigma2[2], FUN = qnorm),
  qtp3(U_sample[,3], mu[3], sigma1[3], sigma2[3], FUN = qnorm),
  qtp3(U_sample[,4], mu[4], sigma1[4], sigma2[4], FUN = qnorm)
)

theta_sample <- exp(eta_sample)
colnames(theta_sample) <- c("lambda","kappa","alpha","beta")
summary(theta_sample)


plot(density(exp(sim[,1])), main = expression(lambda[1]),
     col = "blue", lwd = 2, xlab = expression(lambda[1]), ylab = "Density", ylim = c(0,2))
lines(density(theta_sample[,"lambda"]), col = "red", lwd = 2)
legend("topright", c("MCMC", "Copula-MLE"), col = c("blue","red"), lwd = 2)

plot(density(exp(sim[,2])), main = expression(kappa[1]),
     col = "blue", lwd = 2, xlab = expression(kappa[1]), ylab = "Density")
lines(density(theta_sample[,"kappa"]), col = "red", lwd = 2)
legend("topright", c("MCMC", "Copula-MLE"), col = c("blue","red"), lwd = 2)

plot(density(exp(sim[,3])), main = expression(alpha[1]),
     col = "blue", lwd = 2, xlab = expression(alpha[1]), ylab = "Density", xlim = c(0,10))
lines(density(theta_sample[,"alpha"]), col = "red", lwd = 2)
legend("topright", c("MCMC", "Copula-MLE"), col = c("blue","red"), lwd = 2)

plot(density(exp(sim[,4])), main = expression(beta[1]),
     col = "blue", lwd = 2, xlab = expression(beta[1]), ylab = "Density", xlim = c(-0.5, 10))
lines(density(theta_sample[,"beta"]), col = "red", lwd = 2)
legend("topright", c("MCMC", "Copula-MLE"), col = c("blue","red"), lwd = 2)


