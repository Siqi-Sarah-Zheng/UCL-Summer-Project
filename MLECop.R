rm(list=ls())
library(copula)
library(twopiece)
library(mvtnorm)

sim <- rmvnorm(n = 1000, mean = rep(2,4), sigma = diag(4))


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

