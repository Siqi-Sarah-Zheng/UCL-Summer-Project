##############################
# Normal Approximation
##############################

library(numDeriv)

MAP = nlminb(inits, log_postHR, control = list(iter.max = 1e4))$par


hess = hessian(func = log_postHR, x = MAP)


Sigma = solve(hess)

library(mvtnorm)

set.seed(123)
post_napp = rmvnorm(n = 10000, mean = MAP, sigma = Sigma)

##############################################################
# Graphs for comparisons between mcmc and normal approximation
##############################################################

    
# lambda
plot(density(exp(post_napp[,1])), main= expression(lambda ~ ": MCMC vs Normal Approximation"),
     lwd = 2, col = "blue", xlab=expression(lambda), ylab = "Density")
lines(density(lambdapHR), lwd = 2, col = "red")
legend("topright", c("Normal","MCMC"), col=c("blue","red"), lwd=2)

# kappa
plot(density(exp(post_napp[,2])), main = expression(kappa ~ ": MCMC vs Normal Approximation"),
     lwd = 2, col = "blue", xlab = expression(kappa), ylab = "Density")
lines(density(kappapHR), lwd = 2, col = "red")
legend("topright", c("Normal","MCMC"), col=c("blue","red"), lwd=2)     

# alpha
plot(density(exp(post_napp[,3])), main = expression(alpha ~ ": MCMC vs Normal Approximation"),
     lwd = 2, col = "blue", xlab = expression(alpha), ylab = "Density", xlim = c(0,11))
lines(density(alphapHR), lwd = 2, col = "red")
legend("topright", c("Normal","MCMC"), col=c("blue","red"), lwd=2)  

# beta
plot(density(exp(post_napp[,4])), main = expression(beta ~ ": MCMC vs Normal Approximation"),
     lwd = 2, col = "blue", xlab = expression(beta), ylab = "Density", xlim = c(1,10))
lines(density(betapHR), lwd = 2, col = "red")
legend("topright", c("Normal","MCMC"), col=c("blue","red"), lwd=2)  


