## This script is about unbiasedly estimating the expectation of the limiting
## distribution of an autoregressive process, using coupled chains. We compare
## the "asymptotic inefficiency" defined as variance times expected cost with
## the inefficiency of ergodic averages. The example is taken from "Unbiased
## Monte Carlo: Posterior estimation for intractable/infinite-dimensional
## models", Agapiou, Roberts & Vollmer

library(unbiasedmcmc)
library(lubridate)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)

rm(list = ls())
set.seed(12)

## naive maximal coupling of two Normals
## samples from Normal(mu1, D) and Normal(mu2, D)
maxcoupling <- function(mu1, mu2, sqrtD){
  x <- rnorm(1, mu1, sqrtD)
  if (dnorm(x, mu1, sqrtD, log = TRUE) + log(runif(1)) < dnorm(x, mu2, sqrtD, log = TRUE)){
    return(list(result1 = x, result2 = x, samesame = TRUE))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rnorm(1, mu2, sqrtD)
      reject <- (dnorm(y, mu2, sqrtD, log = TRUE) + log(runif(1)) < dnorm(y, mu1, sqrtD, log = TRUE))
    }
    return(list(result1 = x, result2 = y, samesame = FALSE))
  }
}

## reflection maximal coupling of two Normals 
## samples from Normal(mu1, D) and Normal(mu2, D)
reflmaxcoupling <- function(mu1, mu2, sqrtD, kappa = 1){
  ## first, check if mu1 and mu2 are exactly the same,
  ## because then we can just sample from Normal(mu1, D) and copy the result
  if (identical(mu1, mu2)){
    result1 <- mu1 + sqrtD * rnorm(1, 0, 1)
    result2 <- result1
    samesame <- TRUE
    return(list(result1 = result1, result2 = result2, samesame = samesame))
  }
  dim_ <- length(mu1)
  noise1 <- rnorm(dim_, 0, 1)
  noise2 <- rep(0, dim_)
  logu <- log(runif(1))
  z <- (mu1 - mu2) / sqrtD
  normz <- sqrt(sum(z^2))
  evector <- z / normz
  edotxi <- sum(evector * noise1)
  if (logu < (dnorm(edotxi + kappa * normz, 0, 1, log = TRUE) - dnorm(edotxi, 0, 1, log = TRUE))){
    noise2 <- noise1 + kappa * z
    samesame <- TRUE
  } else {
    noise2 <- noise1 - 2 * edotxi * evector
    samesame <- FALSE
  }
  result1 <- mu1 + noise1 * sqrtD
  result2 <- mu2 + noise2 * sqrtD
  return(list(result1 = result1, result2 = result2, samesame = samesame))
}

## define auto-regressive model
## autoregressive coefficient
rho <- 0.95
## standard deviation
stddev <- sqrt(1-rho^2)
## initial distribution
rinit <- function(){
  return(list(chain_state = 0.))
}
## auto-regressive kernel 
ar_kernel <- function(state){
  return(list(chain_state = rho * state$chain_state + stddev * rnorm(1, mean = 0, sd = 1)))
}

## simulate forward
niterations <- 1e5
state <- rinit()
ar_chain <- rep(0, niterations)
for (iter in 1:niterations){
  state <- ar_kernel(state)
  ar_chain[iter] <- state$chain_state  
}
plot(ar_chain[1:100], type = "l")
hist(ar_chain[200:niterations], prob = TRUE, nclass = 100)
curve(dnorm(x), add = TRUE, lwd = 2, col = "red")

coda::spectrum0.ar(ar_chain[200:niterations])$spec
## this should be approximately equal to 
## sum_t sum_s covar(X_s, X_t)
## 1 + 2 sum_{t\geq 1} covar(X_1, X_t)
## 1 + 2 sum_{t\geq 1} rho^t
## 1 + 2 * rho / (1-rho)
1 + 2 * rho / (1-rho)

## max coupling of AR kernel 
maxcoupled_ar_kernel <- function(state1, state2){
  coupled_states <- maxcoupling(mu1 = rho * state1$chain_state, mu2 = rho * state2$chain_state, sqrtD = stddev)
  return(list(state1 = list(chain_state = coupled_states$result1), 
              state2 = list(chain_state = coupled_states$result2), 
              identical = coupled_states$samesame))
}

## refl-max coupling of AR kernel 
reflmaxcoupled_ar_kernel <- function(state1, state2){
  coupled_states <- reflmaxcoupling(mu1 = rho * state1$chain_state, mu2 = rho * state2$chain_state, sqrtD = stddev)
  return(list(state1 = list(chain_state = coupled_states$result1), 
              state2 = list(chain_state = coupled_states$result2), 
              identical = coupled_states$samesame))
}

state1 = rinit()
state2 = rinit()
maxcoupled_ar_kernel(state1, state2)
## unbiased estimation 
## first, draw meeting times
nrep <- 1e4
maxcoupling_meetingtimes <- foreach(irep = 1:nrep) %dorng% {sample_meetingtime(ar_kernel, maxcoupled_ar_kernel, rinit)$meetingtime}
reflmaxcoupling_meetingtimes <- foreach(irep = 1:nrep) %dorng% {sample_meetingtime(ar_kernel, reflmaxcoupled_ar_kernel, rinit)$meetingtime}
#
summary(unlist(maxcoupling_meetingtimes))
summary(unlist(reflmaxcoupling_meetingtimes))
hist(unlist(maxcoupling_meetingtimes), prob = TRUE, col = rgb(0,0,0), nclass = 100)
hist(unlist(reflmaxcoupling_meetingtimes), prob = TRUE, add = TRUE, nclass = 100, col = rgb(1,0,0,0.5))

## so we can choose e.g. 
lag <- 50
k <- 100
m <- 500
## and produce unbiased estimators as follows
max_uestimators <- foreach(irep = 1:nrep) %dorng% {
  sample_unbiasedestimator(ar_kernel, maxcoupled_ar_kernel, rinit, h = function(x) x, k = k, m = m, lag = lag)
}
reflmax_uestimators <- foreach(irep = 1:nrep) %dorng% {
  sample_unbiasedestimator(ar_kernel, reflmaxcoupled_ar_kernel, rinit, h = function(x) x, k = k, m = m, lag = lag)
}
## check results
mean(sapply(max_uestimators, function(l) l$uestimator))
mean(sapply(reflmax_uestimators, function(l) l$uestimator))

## MSE work product = variance work product
var(sapply(max_uestimators, function(l) l$uestimator)) * mean(sapply(max_uestimators, function(l) l$cost))
var(sapply(reflmax_uestimators, function(l) l$uestimator)) * mean(sapply(reflmax_uestimators, function(l) l$cost))
## compared to 
1 + 2 * rho / (1-rho)
## we lose little in terms of inefficiency if lag, k and m are chosen appropriately

## variance of the unbiased estimator
var(sapply(reflmax_uestimators, function(l) l$uestimator))
## variance of the MCMC estimator with k-1 steps as burn-in and m steps in total 
var(sapply(reflmax_uestimators, function(l) l$mcmcestimator))



