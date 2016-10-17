# In R, the optim function has several options for the optimization method to use, 
# but doesn't offer Newton's method as an option. Newton's method should have a quadratic
# convergence rate and thus be faster than all the optim() methods, so here I write an
# implementation of Newton's method and compare it to the optim() methods by
# doing maximum likelihood estimation on the Weibull distribution.



# This function takes in a sample, y, and parameter values lambda and theta. It outputs
# the Weibull log-likeliood, the gradient, and the Hessian (matrix of second derivatives). 
# The outputs of this function are needed for Newton's Method
weib_loglik <- function(y,lambda,theta) {
    n <- length(y)
    l <- lambda 
    t <- theta
    loglik <- n*log(t) - n*t*log(l) - l^(-t)*sum(y^t) + (t-1)*sum(log(y))
    dl <- -n*t/l + t*l^(-t-1)*sum(y^t)
    dt <- n/t - n*log(l) + sum(log(y)) - sum((y/l)^t*log(y/l))
    dl2 <- n*t/l^2 - t*(t+1)*l^(-t-2)*sum(y^t)
    dt2 <- -n/t^2 - sum((y/l)^t*(log(y/l))^2)
    dldt <- -n/l - l^(-t-1)*sum( y^t * (t*log(l) - t*log(y) - 1) )
    list(Log_Likelihood = loglik, Gradient = c(dl, dt), Hessian = rbind( c(dl2,dldt), c(dldt,dt2) ) )
}


# This function implements Newton's method using the outputs of the above function.
# It stops if the gradient's norm is less than 10^(-8) or if it's done [maxits] iterations


weib_NRmax <- function(y, start_lambda, start_theta, maxits) {
    parms <- c(start_lambda, start_theta)
    P <- cbind(parms)
    for ( i in 1:maxits ) {
        loglik <- weib_loglik(y, parms[1], parms[2])
        g <- loglik$Gradient
        if (as.numeric(sqrt(crossprod(g))) < 10^(-8)) break
        H <- loglik$Hessian
        parms <- parms - solve(H,g)
        P <- cbind(P,as.matrix(parms))
    }
    num_its <- dim(P)[2]
    list(result = P[,num_its], Iterations = num_its)
}

# let's see if it works by generating a sample of 1000 observations from a Weibull(1,1) dist.

y <- rweibull(1000,1,1)

# Try a few different starting values

start1 <- weib_NRmax(y,1,1,100)
start1
start2 <- weib_NRmax(y,1,1/2,100)
start2
start3 <- weib_NRmax(y,1/2,1,100)
start3
start4 <- weib_NRmax(y,0.01,0.01,50)
start4


# Looks like it works 




# Here I run optim() with each method to maximize the likelihood, and make a table (Num_Its)
# of the number of iterations each method required to get the solution.
# As expected, Newton's method gets the solution in fewer iterations than any of the 
# methods available through optim().


# The optim function finds a minimum, so I will tell it to optimize the negative log likelihood

negloglik <- function(start) {
    n <- length(y)
    l <- start[1]
    t <- start[2]
    -(n*log(t) - n*t*log(l) - l^(-t)*sum(y^t) + (t-1)*sum(log(y)))
}

Num_Its <- matrix(0,5,3)
meth <-  c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Newton's")
startval <- cbind(c(1,1),c(1/2,1), c(1,1/2))
y <- rweibull(1000,1,1)

for (i in 1:5) {
    for (j in 1:3) {
        Num_Its[i,j] <- optim(startval[,j], negloglik, method = meth[i])$counts[[1]]
    }
}
Newton_Iterations <- c(weib_NRmax(y,1,1,100)$Iterations, weib_NRmax(y,1/2,1,100)$Iterations, weib_NRmax(y,1,1/2,100)$Iterations)

Num_Its <- rbind(Num_Its, Newton_Iterations)
Num_Its <- as.table(Num_Its)
rownames(Num_Its) <- meth
colnames(Num_Its) <- c('(1, 1)','(1/2, 1)', '(1, 1/2)')

Num_Its <- Num_Its[order(rowMeans(Num_Its)),] # Order by row means

Num_Its

