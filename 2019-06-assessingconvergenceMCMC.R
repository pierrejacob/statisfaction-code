## Pierre E. Jacob -- June 2019

## This script goes through the three examples in the paper:
## "A simulation approach to convergence rates for Markov chain Monte Carlo algorithms"
## by Mary Kathryn Cowles and Jeffrey S. Rosenthal (Statistics and Computing 8.2 (1998): 115-124)
## This is referred to later as the "CR paper".

## In that paper the authors estimate the mixing time of some Gibbs samplers via numerical methods
## that estimate some constants appearing in minorization and drift conditions.

## We revisit the examples with the coupling approach described in
## "Estimating Convergence of Markov chains with L-Lag Couplings." arXiv preprint arXiv:1905.09971 (2019).
## by Niloy Biswas and Pierre E. Jacob.
## For each example, we recode the original Gibbs samplers, and propose simple couplings of them that we can generate from.
## From there we obtain i.i.d. meeting times and can compute upper bounds on the TV from there.

## We start by defining some generic functions that will be used in all examples,
## before defining more specific functions relating to each of the three examples in the CR paper.

# Always start with rm(list = ls()), to annoy R experts
rm(list = ls())
set.seed(21)
library(doParallel) # there will be lots of stuff that can benefit from parallel calculations
registerDoParallel(cores = detectCores()-2) # to run on multiple cores if available but leaving two cores free out of basic courtesy
library(doRNG) # there will be random numbers!
library(dplyr) # some data frame will be manipulated!
library(ggplot2) # there will be plots!

# actually the rest of the code also relies on lme4, tidyr, truncnorm, gridExtra

theme_set(theme_bw()) # a sober theme and various graphical tweaks
theme_update(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
             axis.title.x = element_text(size = 25, margin=margin(20,0,0,0)),
             axis.title.y = element_text(size = 25, angle = 90, margin = margin(0,20,0,0)),
             legend.text = element_text(size = 20), legend.title = element_text(size = 20), title = element_text(size = 30),
             strip.text = element_text(size = 25), strip.background = element_rect(fill="white"), legend.position = "bottom")


# Sample from maximally coupled variables following p and q
# defined through their (log) probability density function dp and dq
# and their random generator rp and rq
rmax_coupling <- function(rp, dp, rq, dq){
  x <- rp(1)
  if (dp(x) + log(runif(1)) < dq(x)){
    return(c(x,x))
  } else {
    reject <- TRUE
    while (reject){
      y <- rq(1)
      reject <- (dq(y) + log(runif(1)) < dp(y))
    }
    return(c(x,y))
  }
}

# Sample a meeting time, taking as arguments:
# - rinit: a function taking no arguments and returning an initial state
# - single_kernel: a function taking a state and return a new state, following the distribution defined by a Markov kernel
# - coupled_kernel: a function taking two states and returning a list containing two new states, following the distribution defined by a joint Markov kernel
# - lag: an integer (default = 1) indicating the lag between the two chains
# - max_iterations: with while loops come great responsabilities...
meetingtimes <- function(rinit, single_kernel, coupled_kernel, lag = 1, max_iterations = Inf){
  # initialize
  state1 <- rinit(); state2 <- rinit()
  # move first chain
  time <- 0
  for (t in 1:lag){
    time <- time + 1
    state1 <- single_kernel(state1)
  }
  # move two chains until meeting (or until max_iterations)
  meetingtime <- Inf
  while (is.infinite(meetingtime) && (time < max_iterations)){
    time <- time + 1
    # use coupled kernel
    coupledstates <- coupled_kernel(state1, state2)
    state1 <- coupledstates$state1
    state2 <- coupledstates$state2
    # check if meeting happens
    allequal <- all(state1 == state2)
    if (allequal) meetingtime <- time
  }
  return(meetingtime)
}

## Compute upper bound on TV between chain at step t and invariant distribution
## based on meeting times obtained with "lag" (produced with above function)
tv_upper_bound <- function(meetingtimes, lag, t){
  return(mean(pmax(0, ceiling((meetingtimes-lag-t) / lag))))
}

## Compute log-density of inverse Gamma
## at x > 0 and with given parameters alpha, beta, given by
##  alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x
dinversegamma <- function(x, alpha, beta){
  return(alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x)
}

## Sample from inverse Gamma with parametrization as above
rinversegamma <- function(n, alpha, beta){
  return(1/rgamma(n = n, shape = alpha, rate = beta))
}

## Sample from maximally coupled inverse Gamma variables
rinversegamma_coupled <- function(alpha1, alpha2, beta1, beta2){
  return(rmax_coupling(function(n) rinversegamma(n, alpha1, beta1),
                function(x) dinversegamma(x, alpha1, beta1),
                function(n) rinversegamma(n, alpha2, beta2),
                function(x) dinversegamma(x, alpha2, beta2)))
}

## Sample from maximally coupled Normal variables (a reflection-maximal coupling might be better but hey).
rnorm_max_coupling <- function(mu1, mu2, sigma1, sigma2){
  return(rmax_coupling(function(n) rnorm(n, mu1, sigma1),
                function(x) dnorm(x, mu1, sigma1, log = TRUE),
                function(n) rnorm(n, mu2, sigma2),
                function(x) dnorm(x, mu2, sigma2, log = TRUE)))
}

## Now, we're in business!

## First example: "a model related to James-Stein estimators"
## data
Y <- c(0.395, 0.375, 0.355, 0.334, 0.313, 0.313, 0.291, 0.269, 0.247, 0.247, 0.224, 0.224,
       0.224, 0.224, 0.224, 0.200, 0.175, 0.148)
ndata <- length(Y)
## The data are Y_i for i in {1,...,K}, with K = 18.
## The model specifies Y_i ~ Normal(theta_i, V), with V known and fixed to 0.00434
V <- 0.00434
## The prior is theta_i ~ Normal(mu, A), with mu ~ flat prior, and A ~ InverseGamma(a,b)
## with density proportional to exp(-b/x) x^{-a-1}; below a = -1, b = 2
a <- -1; b <- 2
## To target the posterior distribution, we consider the following Gibbs sampler:
#- A given rest: IG(a + (K-1)/2, b + (sum_i=1^K (theta_i - theta_bar)^2)/2)
#- mu given rest: Normal(theta_bar, A / K)
#- theta_i given rest: Normal( (mu * V + Y_i * A) / (V + A), A * V / (V + A))
## ... where theta_bar is the average of the theta_i, for i in {1,...,K}

# we store the parameters as (mu, A, theta_1, ..., theta_K) so the parameter space is of dimension 20
# the initialization is as follows (the initialization of mu and A is irrelevant)
rinit <- function() c(0, 1, rep(mean(Y), ndata))

# sample from Markov kernel
single_kernel <- function(state){
  theta <- state[3:(ndata+2)]
  theta_bar <- mean(state[3:(ndata+2)])
  # update of A given rest
  A <- rinversegamma(1, a + 0.5 * (ndata-1), b + 0.5 * sum((theta - theta_bar)^2))
  # update of mu given rest
  mu <- rnorm(1, theta_bar, sqrt(A/ndata))
  # update of each theta_i
  theta <- rnorm(ndata, (mu * V + Y * A) / (V + A), sqrt(A * V / (V + A)))
  return(c(mu, A, theta))
}

# sample from coupled kernel
coupled_kernel <- function(state1, state2){
  theta1 <- state1[3:(ndata+2)];  theta_bar1 <- mean(theta1)
  theta2 <- state2[3:(ndata+2)];  theta_bar2 <- mean(theta2)
  # update of A given rest
  As <- rinversegamma_coupled(alpha1 = a + 0.5 * (ndata-1), alpha2 = a + 0.5 * (ndata-1),
                              beta1 = b + 0.5 * sum((theta1 - theta_bar1)^2),
                              beta2 = b + 0.5 * sum((theta2 - theta_bar2)^2))
  A1 <- As[1];  A2 <- As[2]
  # update of mu given rest
  mus <- rnorm_max_coupling(theta_bar1, theta_bar2, sqrt(A1/ndata), sqrt(A2/ndata))
  mu1 <- mus[1];  mu2 <- mus[2]
  # # update of each theta_i
  thetas <- matrix(nrow = ndata, ncol = 2)
  for (idata in 1:ndata){
    thetas[idata,] <- rnorm_max_coupling((mu1 * V + Y[idata] * A1) / (V + A1), (mu2 * V + Y[idata] * A2) / (V + A2),
                                         sqrt(A1 * V / (V + A1)), sqrt(A2 * V / (V + A2)))
  }
  theta1 <- thetas[,1]; theta2 <- thetas[,2]
  return(list(state1 = c(mu1, A1, theta1), state2 = c(mu2, A2, theta2)))
}

## Now, how many steps do we need to be close to stationarity in TV?
## Let's find out using meeting times.
NREP <- 1e3
lag_example1 <- 10
meetingtimes_example1 <- foreach(irep = 1:NREP) %dorng% {meetingtimes(rinit, single_kernel, coupled_kernel, lag_example1) }

niterations_example1 <- 6
ubounds_example1 <- sapply(1:niterations_example1, function(t) tv_upper_bound(unlist(meetingtimes_example1), lag_example1, t))
g_tvbounds_examples1 <- qplot(x = 1:niterations_example1, y = ubounds_example1, geom = "line") + scale_x_continuous(breaks = 1:niterations_example1) + ylab("TV upper bounds") + xlab("iteration")
g_tvbounds_examples1 <- g_tvbounds_examples1 + labs(title = "Example 1: baseball") + ylim(0,1)
g_tvbounds_examples1
## it seems that the Gibbs sampler has essentially converged in less than 4 steps
## we can check: let's run many chains for e.g. 4 steps, and compare the posterior marginals
## with approximation obtained with one long chain

## lots of short chains
NREP <- 1e4
nstep <- 4
shortchains <- foreach(irep = 1:NREP, .combine = rbind) %dorng% {
  state <- rinit(); for (istep in 1:nstep)  state <- single_kernel(state)
  state
}

### run one long chain
niterations <- 1e5 # crazy long chain!
burnin <- 200 # ultra conservative burn-in!
longchain <- matrix(nrow = niterations+burnin, ncol = 20)
longchain[1,] <- rinit()
for (iteration in 2:(niterations+burnin)){
  longchain[iteration,] <- single_kernel(longchain[iteration-1,])
}
longchain <- longchain[(burnin+1):(niterations+burnin),]

# compare marginals
g1 <- qplot(x = longchain[,1], geom = "blank") + geom_density() + geom_density(aes(x = shortchains[,1]), linetype = 2) + xlab(expression(mu))
g2 <- qplot(x = longchain[,2], geom = "blank") + geom_density() + geom_density(aes(x = shortchains[,2]), linetype = 2) + xlab(expression(A))
g3 <- qplot(x = longchain[,3], geom = "blank") + geom_density() + geom_density(aes(x = shortchains[,3]), linetype = 2) + xlab(expression(theta[1]))
gridExtra::grid.arrange(g1, g2, g3, nrow = 1)
## indeed the marginal distributions look very similar!

## Second example: "variance components model"
## the Dyestuff data set is available in the package lme4
library(lme4)
# ?Dyestuff
mu_0 <- 0; sigma_0 <- 10^6; sigma_02 <- sigma_0^2
a1 <- 0.5; b1 <- 1
a2 <- 0; b2 <- 0
K <- 6; J <- 5
Ydataframe <- Dyestuff %>% group_by(Batch) %>% mutate(individual = 1:n()) %>% ungroup()
# let's look at the data
Ydataframe %>% group_by(Batch) %>% summarise(mean_per_batch = mean(Yield), var_per_batch = var(Yield))
# manipulate the data to ease forthcoming calculations
Ywide <- tidyr::spread(Ydataframe, Batch, Yield) %>% select(-individual)
Ywide <- t(as.matrix(Ywide)) # data as matrix, with each row corresponding to a batch
Ymean_per_batch <- rowMeans(Ywide)
Ymean <- mean(Ywide)
# note: here we're doing some really advanced stuff: removing a vector of length K to a matrix with K rows
# does a row-wise subtraction with the elements of the vector
# e.g. Ywide - c(1,1,1,2,2,2) removes 1 to the first 3 rows, and removes 2 to the next 3 rows

## Gibbs sampler, as in the CR paper
single_kernel <- function(state){
  # extract parameters from vector state
  sigma_theta2 <- state[1];  sigma_e2 <- state[2];  mu <- state[3];  theta <- state[4:(4+K-1)]
  # update of sigma_theta2
  sigma_theta2 <- rinversegamma(1, a1 + 0.5 * K, b1 + 0.5 * sum((theta - mu)^2))
  # update of sigma_22
  sigma_e2 <- rinversegamma(1, a2 + 0.5 * K * J, b2 + 0.5 * sum((Ywide - theta)^2))
  # update of mu
  mean_mu <- (sigma_theta2 * mu_0 + sigma_02 * sum(theta)) / (sigma_theta2 + K * sigma_02)
  var_mu <- (sigma_theta2*sigma_02)/(sigma_theta2 + K * sigma_02)
  mu <- rnorm(1, mean = mean_mu, sd = sqrt(var_mu))
  # update of each theta
  mean_theta <- (J * sigma_theta2 * Ymean_per_batch + sigma_e2 * mu) / (J * sigma_theta2 + sigma_e2)
  var_theta <- (sigma_theta2 * sigma_e2) / (J * sigma_theta2 + sigma_e2)
  theta <- rnorm(K, mean = mean_theta, sd = sqrt(var_theta))
  return(c(sigma_theta2, sigma_e2, mu, theta))
}
## Coupled Gibbs sampler, with maximally coupled updates
coupled_kernel <- function(state1, state2){
  sigma_theta21 <- state1[1];  sigma_e21 <- state1[2];  mu1 <- state1[3];  theta1 <- state1[4:(4+K-1)]
  sigma_theta22 <- state2[1];  sigma_e22 <- state2[2];  mu2 <- state2[3];  theta2 <- state2[4:(4+K-1)]
  # update of sigma_theta2
  sigma_theta2_ <- rinversegamma_coupled(a1 + 0.5 * K, a1 + 0.5 * K,
                                         b1 + 0.5 * sum((theta1 - mu1)^2),
                                         b1 + 0.5 * sum((theta2 - mu2)^2))
  sigma_theta21 <- sigma_theta2_[1]; sigma_theta22 <- sigma_theta2_[2]
  # update of sigma_22
  sigma_e2_ <- rinversegamma_coupled(a2 + 0.5 * K * J, a2 + 0.5 * K * J,
                                     b2 + 0.5 * sum((Ywide - theta1)^2),
                                     b2 + 0.5 * sum((Ywide - theta2)^2))
  sigma_e21 <- sigma_e2_[1]; sigma_e22 <- sigma_e2_[2]
  # update of mu
  mean_mu1 <- (sigma_theta21 * mu_0 + sigma_02 * sum(theta1)) / (sigma_theta21 + K * sigma_02)
  var_mu1  <- (sigma_theta21*sigma_02)/(sigma_theta21 + K * sigma_02)
  mean_mu2 <- (sigma_theta22 * mu_0 + sigma_02 * sum(theta2)) / (sigma_theta22 + K * sigma_02)
  var_mu2  <- (sigma_theta22*sigma_02)/(sigma_theta22 + K * sigma_02)
  mu_ <- rnorm_max_coupling(mean_mu1, mean_mu2, sqrt(var_mu1), sqrt(var_mu2))
  mu1 <- mu_[1]; mu2 <- mu_[2]
  # update of each theta
  mean_theta1 <- (J * sigma_theta21 * Ymean_per_batch + sigma_e21 * mu1) / (J * sigma_theta21 + sigma_e21)
  var_theta1 <- (sigma_theta21 * sigma_e21) / (J * sigma_theta21 + sigma_e21)
  mean_theta2 <- (J * sigma_theta22 * Ymean_per_batch + sigma_e22 * mu2) / (J * sigma_theta22 + sigma_e22)
  var_theta2 <- (sigma_theta22 * sigma_e22) / (J * sigma_theta22 + sigma_e22)
  for (k in 1:K){
    theta_k <- rnorm_max_coupling(mean_theta1[k], mean_theta2[k], sqrt(var_theta1), sqrt(var_theta2))
    theta1[k] <- theta_k[1]; theta2[k] <- theta_k[2]
  }
  return(list(state1 = c(sigma_theta21, sigma_e21, mu1, theta1),
              state2 = c(sigma_theta22, sigma_e22, mu2, theta2)))
}

## We can re-implement the V function of the CR paper
v1 <- mean((Ywide - Ymean_per_batch)^2)
v2 <- mean((Ymean_per_batch - Ymean)^2)
rinit <- function() c(1, 1, Ymean, (J * v1 * Ymean_per_batch + v2 * Ymean) / (J * v1 + v2))
V_function <- function(state){
  sigma_theta2 <- state[1];  sigma_e2 <- state[2];  mu <- state[3];  theta <- state[4:(4+K-1)]
  return(mean((theta - (J * v1 * Ymean_per_batch + v2 * Ymean) / (J * v1 + v2))^2) + (mu  - Ymean)^2)
}
V_function(rinit()) == 0
## So, it seems we got the same initial distribution as in the CR paper

## Now, how many steps do we need to be close to stationarity in TV??
NREP <- 1e3
lag_example2 <- 500
meetingtimes_example2 <- foreach(irep = 1:NREP) %dorng% { meetingtimes(rinit, single_kernel, coupled_kernel, lag_example2)}

niterations_example2 <- 5e2
ubounds_example2 <- sapply(1:niterations_example2, function(t) tv_upper_bound(unlist(meetingtimes_example2), lag_example2, t))
g_tvbounds_examples2 <- qplot(x = 1:niterations_example2, y = ubounds_example2, geom = "line")
g_tvbounds_examples2 <- g_tvbounds_examples2 + ylab("TV upper bounds") + xlab("iteration") + labs(title = "Example 2: dyestuff") + ylim(0,1)
g_tvbounds_examples2


## it seems that the Gibbs sampler has essentially converged in less than 1000 steps

## we can check: let's run many chains for e.g. 500 steps, and compare the posterior marginals
## with an approximation obtained with one long chain
## lots of short chains
NREP <- 1e4
nstep <- 500
shortchains <- foreach(irep = 1:NREP, .combine = rbind) %dorng% {
  state <- rinit()
  for (istep in 1:nstep)  state <- single_kernel(state)
  state
}

### one long chain
burnin <- 1e4 # conservative burn-in
niterations <- 1e5 # long chain
longchain <- matrix(nrow = niterations+burnin, ncol = 9)
longchain[1,] <- rinit()
for (iteration in 2:(niterations+burnin)){
  longchain[iteration,] <- single_kernel(longchain[iteration-1,])
}
longchain <- longchain[(burnin+1):(niterations+burnin),]

#
g1 <- qplot(x = longchain[,1], geom = "blank") + geom_density() + geom_density(aes(x = shortchains[,1]), linetype = 2) + xlab(expression(sigma[theta]^2)) + scale_x_log10()
g2 <- qplot(x = longchain[,2], geom = "blank") + geom_density() + geom_density(aes(x = shortchains[,2]), linetype = 2) + xlab(expression(sigma[e]^2)) + scale_x_log10()
g3 <- qplot(x = longchain[,4], geom = "blank") + geom_density() + geom_density(aes(x = shortchains[,4]), linetype = 2) + xlab(expression(theta[1]))
gridExtra::grid.arrange(g1, g2, g3, nrow = 1)
## indeed the marginal distributions look very similar!!

## Third example: "ordinal probit model"
## This example stems from "Accelerating Monte Carlo Markov chain convergence for cumulative-link generalized linear models" by Mary Kathryn Cowles
## It is based on simulated data so might not be exactly as in the CR paper or the above reference. But hopefully it's close enough?

## generate data with n=50
nobs <- 50
beta0_dgp <- +1
beta1_dgp <- -2
X <- rnorm(nobs) # Normal covariates
ystar_dgp <- beta0_dgp + beta1_dgp * X + rnorm(nobs) # generate ystar
w <- as.numeric(cut(ystar_dgp, 3, labels = c(-1,0,+1)))-2 # generate w based on thirds of ystar
X <- cbind(1, X) # add intercept
XtXinverse <- solve(t(X) %*% X) # pre-computation of (X'X)^{-1}

# initial distribution; didn't find the specification in the CR paper, so I made one up
rinit <- function(){
  # we need ystar to satisfy ystar[w==-1] < ystar[w==0] < ystar[w==1]
  ystar <- rep(0, nobs)
  ystar[w==-1] <- runif(sum(w==-1), min = -3/2, max = -1/2)
  ystar[w== 0] <- runif(sum(w== 0), min = -1/2, max = +1/2)
  ystar[w==+1] <- runif(sum(w==+1), min = +1/2, max = +3/2)
  return(c(rnorm(2), runif(1, min = 0.05, max = 10), ystar))
}

# Gibbs sampler as in the CR paper
single_kernel <- function(state){
  beta0 <- state[1]; beta1 <- state[2]; cutoff <- state[3]; ystar <- state[4:(4+nobs-1)]
  B <- XtXinverse %*% (t(X) %*% ystar)
  # update beta
  betas <- mvtnorm::rmvnorm(1, B, XtXinverse)
  beta0 <- betas[1]; beta1 <- betas[2]
  # update cutoff
  M0 <- max(ystar[w == 0])
  m1 <- min(ystar[w == 1])
  cutoff <- runif(1, max(M0, 0.05), m1)
  # update ystar, using truncnorm package
  ystar[w==-1] <- truncnorm::rtruncnorm(n = sum(w==-1), mean = beta0 + beta1 * X[w==-1,2], sd = 1, a = -Inf, b = 0)
  ystar[w== 0] <- truncnorm::rtruncnorm(n = sum(w== 0), mean = beta0 + beta1 * X[w== 0,2], sd = 1, a = 0, b = cutoff)
  ystar[w==+1] <- truncnorm::rtruncnorm(n = sum(w==+1), mean = beta0 + beta1 * X[w==+1,2], sd = 1, a = cutoff, b = +Inf)
  return(c(beta0, beta1, cutoff, ystar))
}

### function to sample from max coupling of two multivariate Normals with same variance
rmvnorm_coupling <- function(mu1, mu2, Sigma){
  x <- mvtnorm::rmvnorm(1, mu1, Sigma)
  if (mvtnorm::dmvnorm(x, mu1, Sigma, log = TRUE) + log(runif(1)) < mvtnorm::dmvnorm(x, mu2, Sigma, log = TRUE)){
    return(rbind(x,x))
  } else {
    reject <- TRUE
    while (reject){
      y <- mvtnorm::rmvnorm(1, mu2, Sigma)
      reject <- (mvtnorm::dmvnorm(y, mu2, Sigma, log = TRUE) + log(runif(1)) < mvtnorm::dmvnorm(y, mu1, Sigma, log = TRUE))
    }
    return(rbind(x,y))
  }
}

## coupled Gibbs kernel
coupled_kernel <- function(state1, state2){
  beta01 <- state1[1]; beta11 <- state1[2]; cutoff1 <- state1[3]; ystar1 <- state1[4:(4+nobs-1)]
  beta02 <- state2[1]; beta12 <- state2[2]; cutoff2 <- state2[3]; ystar2 <- state2[4:(4+nobs-1)]
  B1 <- XtXinverse %*% (t(X) %*% ystar1)
  B2 <- XtXinverse %*% (t(X) %*% ystar2)
  # update beta
  betascoupled <- rmvnorm_coupling(B1, B2, XtXinverse)
  betas1 <- betascoupled[1,]; betas2 <- betascoupled[2,]
  beta01 <- betas1[1]; beta11 <- betas1[2]
  beta02 <- betas2[1]; beta12 <- betas2[2]
  # update cutoff
  M01 <- max(ystar1[w == 0])
  m11 <- min(ystar1[w == 1])
  M02 <- max(ystar2[w == 0])
  m12 <- min(ystar2[w == 1])
  cutoffs_ <- rmax_coupling(rp = function(n) runif(n, min = max(M01, 0.05), max = m11),
                            dp = function(x){ dunif(x, min = max(M01, 0.05), max = m11, log = TRUE)},
                            rq = function(n) runif(n, min = max(M02, 0.05), max = m12),
                            dq = function(x){ dunif(x, min = max(M02, 0.05), max = m12, log = TRUE)})
  cutoff1 <- cutoffs_[1]; cutoff2 <- cutoffs_[2]
  # update ystar
  for (iobs in 1:nobs){
    mean1 <- beta01 + beta11 * X[iobs,2]
    mean2 <- beta02 + beta12 * X[iobs,2]
    if (w[iobs] == -1){
      a1 <- -Inf; b1 <- 0
      a2 <- -Inf; b2 <- 0
    }
    if (w[iobs] == 0){
      a1 <- 0; b1 <- cutoff1
      a2 <- 0; b2 <- cutoff2
    }
    if (w[iobs] == 1){
      a1 <- cutoff1; b1 <- +Inf
      a2 <- cutoff2; b2 <- +Inf
    }
    ystar_coupled_iobs <- rmax_coupling(rp = function(n) truncnorm::rtruncnorm(n, mean = mean1, sd = 1, a = a1, b = b1),
                                        dp = function(x) log(truncnorm::dtruncnorm(x, mean = mean1, sd = 1, a = a1, b = b1)),
                                        rq = function(n) truncnorm::rtruncnorm(n, mean = mean2, sd = 1, a = a2, b = b2),
                                        dq = function(x) log(truncnorm::dtruncnorm(x, mean = mean2, sd = 1, a = a2, b = b2)))
    ystar1[iobs] <- ystar_coupled_iobs[1]
    ystar2[iobs] <- ystar_coupled_iobs[2]
  }
  return(list(state1 = c(beta01, beta11, cutoff1, ystar1), state2 = c(beta02, beta12, cutoff2, ystar2)))
}

## Now, how many steps do we need to be close to stationarity in TV??
NREP <- 1e3
lag_example3 <- 500
meetingtimes_example3 <- foreach(irep = 1:NREP) %dorng% {
  meetingtimes(rinit, single_kernel, coupled_kernel, lag_example3)
}

niterations_example3 <- 4e2
ubounds_example3 <- sapply(1:niterations_example3, function(t) tv_upper_bound(unlist(meetingtimes_example3), lag_example3, t))
g_tvbounds_examples3 <- qplot(x = 1:niterations_example3, y = ubounds_example3, geom = "line")
g_tvbounds_examples3 <- g_tvbounds_examples3 + ylab("TV upper bounds") + xlab("iteration") + labs(title = "Example 3: ordinal probit") + ylim(0,1)
g_tvbounds_examples3


## we can check: let's run many chains for e.g. 300 steps, and compare the posterior marginals
## with an approximation obtained with one long chain
## lots of short chains
NREP <- 1e3
nstep <- 300
shortchains <- foreach(irep = 1:NREP, .combine = rbind) %dorng% {
  state <- rinit()
  for (istep in 1:nstep)  state <- single_kernel(state)
  state
}

### one long chain
burnin <- 1e4
niterations <- 1e5
longchain <- matrix(nrow = niterations+burnin, ncol = nobs+3)
longchain[1,] <- rinit()
for (iteration in 2:(niterations+burnin)){
  longchain[iteration,] <- single_kernel(longchain[iteration-1,])
}
longchain <- longchain[(burnin+1):(niterations+burnin),]

#
g1 <- qplot(x = longchain[,1], geom = "blank") + geom_density() + geom_density(aes(x = shortchains[,1]), linetype = 2) + xlab(expression(beta[0]))
g2 <- qplot(x = longchain[,2], geom = "blank") + geom_density() + geom_density(aes(x = shortchains[,2]), linetype = 2) + xlab(expression(beta[1]))
g3 <- qplot(x = longchain[,3], geom = "blank") + geom_density() + geom_density(aes(x = shortchains[,3]), linetype = 2) + xlab("cut off")
gridExtra::grid.arrange(g1, g2, g3, nrow = 1)
## indeed the marginal distributions look very similar!!


# ggsave(filename = "2019-06-assessing-3.png", plot = gridExtra::grid.arrange(g_tvbounds_examples3, nrow = 1),
#        width = 7, height = 7)
# ggsave(filename = "2019-06-assessing-12.png", plot = gridExtra::grid.arrange(g_tvbounds_examples1, g_tvbounds_examples2, nrow = 1),
#        width = 12, height = 7)

