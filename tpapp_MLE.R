
sim <- chainHR[ind,]

nloglik <- function(par){
  
  mu <- par[1:4]
  sigma1 <- exp(par[5:8])
  sigma2 <- exp(par[9:12])
  
  # Calculate log_q using the pre-calculated grid.
  
  log_q_matrix <- matrix(NA, nrow = nrow(sim), ncol = ncol(sim))
  for(j in 1:nrow(sim)){
  for (i in 1:4) {
    log_q_matrix[j,i] <- dtp3(sim[j,i], mu[i], sigma1[i], sigma2[i], FUN = dnorm, log = TRUE )
  }
  }
  
  out <- -sum(log_q_matrix)
  
return(out)  
  
}


OPT <- nlminb(start = c(rep(0,4),rep(1,8)), objective = nloglik, control = list(iter.max = 10000))


tempf1 <- Vectorize(function(x) dtp3(x, OPT$par[1], exp(OPT$par[5]), exp(OPT$par[9]), FUN = dnorm, log = FALSE ))
tempf2 <- Vectorize(function(x) dtp3(x, OPT$par[2], exp(OPT$par[6]), exp(OPT$par[10]), FUN = dnorm, log = FALSE ))
tempf3 <- Vectorize(function(x) dtp3(x, OPT$par[3], exp(OPT$par[7]), exp(OPT$par[11]), FUN = dnorm, log = FALSE ))
tempf4 <- Vectorize(function(x) dtp3(x, OPT$par[4], exp(OPT$par[8]), exp(OPT$par[12]), FUN = dnorm, log = FALSE ))


plot(density(sim[,1]), main= expression(lambda ~ ": MCMC vs TPNormal Approximation"),
     lwd = 2, col = "blue", xlab=expression(lambda), ylab = "Density")
curve(tempf1,0.1,1.1, col = "red", add= T, lwd = 2, n = 1000)


plot(density(sim[,2]), main= expression(lambda ~ ": MCMC vs TPNormal Approximation"),
     lwd = 2, col = "blue", xlab=expression(lambda), ylab = "Density")
curve(tempf2,-3,-1, col = "red", add= T, lwd = 2, n = 1000)


plot(density(sim[,3]), main= expression(lambda ~ ": MCMC vs TPNormal Approximation"),
     lwd = 2, col = "blue", xlab=expression(lambda), ylab = "Density")
curve(tempf3,0,2.5, col = "red", add= T, lwd = 2, n = 1000)


plot(density(sim[,4]), main= expression(lambda ~ ": MCMC vs TPNormal Approximation"),
     lwd = 2, col = "blue", xlab=expression(lambda), ylab = "Density")
curve(tempf4,-0.5,2.5, col = "red", add= T, lwd = 2, n = 1000)


