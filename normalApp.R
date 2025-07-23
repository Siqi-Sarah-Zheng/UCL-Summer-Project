##############################
# Normal Approximation
##############################

library(numDeriv)

OPTMAP = nlminb(inits, log_postHR, control = list(iter.max = 1e4))

MAP = OPTMAP$par


hess = hessian(func = log_postHR, x = MAP)


Sigma = solve(hess)

library(mvtnorm)

set.seed(123)
post_napp = rmvnorm(n = 1000, mean = MAP, sigma = Sigma)


plot(density(exp(post_napp[,1])))


lnormc <- -OPTMAP$objective + 0.5*length(MAP)*log(2*pi) - 0.5*determinant(hess, logarithm = TRUE)$modulus
lnormc <- as.numeric(lnormc)

