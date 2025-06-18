##############################
# Normal Approximation
##############################

library(numDeriv)

MAP = nlminb(inits, log_postHR, control = list(iter.max = 1e4))$par


hess = hessian(func = log_postHR, x = MAP)


Sigma = solve(hess)

library(mvtnorm)

set.seed(123)
post_napp = rmvnorm(n = 1000, mean = MAP, sigma = Sigma)


plot(density(exp(post_napp[,1])))


