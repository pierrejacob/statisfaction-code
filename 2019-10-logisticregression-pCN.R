## This script is about unbiased MCMC, with a pCN proposal in a
## Metropolis--Hastings algorithm targeting the posterior distribution in a
## simulated logistic regression. A standard MH algorithm is also implemented
## for comparison The example is taken from "Unbiased Monte Carlo: Posterior
## estimation for intractable/infinite-dimensional models", Agapiou, Roberts &
## Vollmer

library(unbiasedmcmc)
library(lubridate)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
rm(list = ls())
set.seed(12)

### logistic regression setting 
target_dimension <- 3
N <- 100
X <- cbind(matrix(rnorm((target_dimension-1) * N), nrow = N, ncol = target_dimension-1), rep(1, N))
## generate the Y
beta_dgp <- c(1,2,3)
##
Y <- runif(N) < (1 / (1 + exp(- X %*% beta_dgp)[,1]))

library(Rcpp)
## function to compute log-likelihood in logistic regression
cppFunction("
double loglikelihood_(const NumericVector & Y,
                      const NumericMatrix & X, const NumericVector & theta){
  double logpr = 0.;
  int p = theta.size();
  int ny = X.rows();
  for (int i = 0; i < ny; i++){
    double thetaX = 0;
    for (int j = 0; j < p; j++){
      thetaX = thetaX + theta(j) * X(i,j);
    }
    logpr = logpr + (Y(i)-1) * thetaX - log(1+exp(-thetaX));
  }
  return(logpr);
}")

## define log pdf of posterior distribution
## assumes a Normal(0,I) prior on parameters
logtarget <- function(chain_state) loglikelihood_(Y, X, chain_state) + sum(dnorm(chain_state, log = TRUE))
## find MAP estimate
map <- optim(par = c(0,1,2), fn = function(v) - logtarget(v))$par
## and Hessian at the MAP estimate
precision <- numDeriv::hessian(function(v) - logtarget(v), map)
c <- map
C <- solve(precision)
## The posterior might be close to Normal(c,C)
## 
Sigma_proposal <- C
Sigma_proposal_chol <- chol(Sigma_proposal)
Sigma_proposal_chol_inv <- solve(chol(Sigma_proposal))
zeromean <- rep(0, target_dimension)

## Initial distribution Normal(c,C)
rinit <- function(){
  chain_state <- fast_rmvnorm_chol(1, mean = c, chol = Sigma_proposal_chol)[1,]
  current_pdf <- logtarget(chain_state)
  return(list(chain_state = chain_state, current_pdf = current_pdf))
}

## MH kernel with random walk proposals
rwmh_kernel <- function(state){
  chain_state <- state$chain_state
  current_pdf <- state$current_pdf
  proposal_value <- chain_state + fast_rmvnorm_chol(1, zeromean, Sigma_proposal_chol)
  proposal_pdf <- logtarget(proposal_value)
  accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
  if (accept){
    return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
  } else {
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
}

## Coupled kernel 
coupled_rwmh_kernel <- function(state1, state2){
  chain_state1 <- state1$chain_state; current_pdf1 <- state1$current_pdf
  chain_state2 <- state2$chain_state; current_pdf2 <- state2$current_pdf
  proposal_value <- rmvnorm_reflectionmax(chain_state1, chain_state2, Sigma_proposal_chol, Sigma_proposal_chol_inv)
  proposal1 <- proposal_value$xy[,1]; proposal_pdf1 <- logtarget(proposal1)
  if (proposal_value$identical){
    proposal2 <- proposal1; proposal_pdf2 <- proposal_pdf1
  } else {
    proposal2 <- proposal_value$xy[,2]; proposal_pdf2 <- logtarget(proposal2)
  }
  logu <- log(runif(1))
  accept1 <- FALSE; accept2 <- FALSE
  if (is.finite(proposal_pdf1)){
    accept1 <- (logu < (proposal_pdf1 - current_pdf1))
  }
  if (is.finite(proposal_pdf2)){
    accept2 <- (logu < (proposal_pdf2 - current_pdf2))
  }
  if (accept1){
    chain_state1 <- proposal1
    current_pdf1 <- proposal_pdf1
  }
  if (accept2){
    chain_state2 <- proposal2
    current_pdf2 <- proposal_pdf2
  }
  identical_ <- proposal_value$identical && accept1 && accept2
  return(list(state1 = list(chain_state = chain_state1, current_pdf = current_pdf1),
              state2 = list(chain_state = chain_state2, current_pdf = current_pdf2),
              identical = identical_))
}


state <- rinit()
niterations <- 1e6
rwmh_chain <- matrix(nrow = niterations, ncol = target_dimension)
for (iter in 1:niterations){
  state <- rwmh_kernel(state)
  rwmh_chain[iter,] <- state$chain_state
}

## now MH with pCN proposals
rho <- 0.8
sqrt1mrho2 <- sqrt(1 - rho^2)
dimension <- dim(Sigma_proposal)[1]
# single kernel
pcnmh_kernel <- function(state){
  chain_state <- matrix(state$chain_state, nrow = 1)
  current_pdf <- state$current_pdf
  proposal_value <- fast_rmvnorm_chol(1, c + rho * (chain_state-c), sqrt1mrho2 * Sigma_proposal_chol)
  proposal_pdf <- logtarget(proposal_value)
  acceptratio <- (proposal_pdf - current_pdf) + 
    fast_dmvnorm_chol_inverse(chain_state, c + rho * (proposal_value-c), Sigma_proposal_chol_inv / sqrt1mrho2) - 
    fast_dmvnorm_chol_inverse(proposal_value, c + rho * (chain_state-c), Sigma_proposal_chol_inv / sqrt1mrho2)
  accept <- (log(runif(1)) < acceptratio)
  if (accept){
    return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
  } else {
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
}

pcnmh_chain <- matrix(nrow = niterations, ncol = target_dimension)
state <- rinit()
for (iter in 1:niterations){
  state <- pcnmh_kernel(state)
  pcnmh_chain[iter,] <- state$chain_state
}

## compare posterior mean approximations
colMeans(rwmh_chain[1000:niterations,])
colMeans(pcnmh_chain[1000:niterations,])

## compare posterior variance approximations
cov(rwmh_chain[1000:niterations,])
cov(pcnmh_chain[1000:niterations,])

par(mfrow = c(1,3))
hist(rwmh_chain[1000:niterations,1], prob = TRUE, nclass = 100, col = rgb(0,0,0), main = "component 1", xlab = "x1")
hist(pcnmh_chain[1000:niterations,1], prob = TRUE, nclass = 100, col = rgb(1,0,0,0.5), add = TRUE)
abline(v = map[1])
       
hist(rwmh_chain[1000:niterations,2], prob = TRUE, nclass = 100, col = rgb(0,0,0), main = "component 2", xlab = "x2")
hist(pcnmh_chain[1000:niterations,2], prob = TRUE, nclass = 100, col = rgb(1,0,0,0.5), add = TRUE)
abline(v = map[2])

hist(rwmh_chain[1000:niterations,3], prob = TRUE, nclass = 100, col = rgb(0,0,0), main = "component 3", xlab = "x3")
hist(pcnmh_chain[1000:niterations,3], prob = TRUE, nclass = 100, col = rgb(1,0,0,0.5), add = TRUE)
abline(v = map[3])

# coupled kernel
coupled_pcnmh_kernel <- function(state1, state2){
  chain_state1 <- matrix(state1$chain_state, nrow = 1); current_pdf1 <- state1$current_pdf
  chain_state2 <- matrix(state2$chain_state, nrow = 1); current_pdf2 <- state2$current_pdf
  proposal_value <- rmvnorm_reflectionmax(c + rho * (chain_state1-c), c + rho * (chain_state2-c), sqrt1mrho2 * Sigma_proposal_chol, Sigma_proposal_chol_inv / sqrt1mrho2)
  proposal1 <- t(proposal_value$xy[,1,drop=F]); proposal_pdf1 <- logtarget(proposal1)
  if (proposal_value$identical){
    proposal2 <- proposal1; proposal_pdf2 <- proposal_pdf1
  } else {
    proposal2 <- t(proposal_value$xy[,2,drop=F]); proposal_pdf2 <- logtarget(proposal2)
  }
  logu <- log(runif(1))
  accept1 <- FALSE; accept2 <- FALSE
  if (is.finite(proposal_pdf1)){
    cat(proposal_pdf1, "vs", current_pdf1, "\n")
    acceptratio1 <- (proposal_pdf1 - current_pdf1) + 
      fast_dmvnorm_chol_inverse(chain_state1, c + rho * (proposal1-c), Sigma_proposal_chol_inv / sqrt1mrho2) - 
      fast_dmvnorm_chol_inverse(proposal1, c + rho * (chain_state1-c), Sigma_proposal_chol_inv / sqrt1mrho2)
    accept1 <- (logu < acceptratio1)
  }
  if (is.finite(proposal_pdf2)){
    acceptratio2 <- (proposal_pdf2 - current_pdf2) + 
      fast_dmvnorm_chol_inverse(chain_state2, c + rho * (proposal2-c), Sigma_proposal_chol_inv / sqrt1mrho2) - 
      fast_dmvnorm_chol_inverse(proposal2, c + rho * (chain_state2-c), Sigma_proposal_chol_inv / sqrt1mrho2)
    accept2 <- (logu < acceptratio2)
  }
  if (accept1){
    chain_state1 <- proposal1
    current_pdf1 <- proposal_pdf1
  }
  if (accept2){
    chain_state2 <- proposal2
    current_pdf2 <- proposal_pdf2
  }
  identical_ <- proposal_value$identical && accept1 && accept2
  return(list(state1 = list(chain_state = chain_state1, current_pdf = current_pdf1),
              state2 = list(chain_state = chain_state2, current_pdf = current_pdf2),
              identical = identical_))
}

## generate meeting times
nrep <- 1e4
rwmh_meetingtimes <- foreach(irep = 1:nrep) %dorng% {sample_meetingtime(rwmh_kernel, coupled_rwmh_kernel, rinit)$meetingtime}
pcnmh_meetingtimes <- foreach(irep = 1:nrep) %dorng% {sample_meetingtime(pcnmh_kernel, coupled_pcnmh_kernel, rinit)$meetingtime}
##
summary(unlist(rwmh_meetingtimes))
summary(unlist(pcnmh_meetingtimes))
## so we can choose e.g. lag = 50, k = 100, m = 500
lag <- 50
k <- 100
m <- 500
## generate unbiased estimators
rwmh_uestimators <- foreach(irep = 1:nrep) %dorng% sample_unbiasedestimator(rwmh_kernel, coupled_rwmh_kernel, rinit, h = function(x) x, k = k, m = m, lag = lag)
pcnmh_uestimators <- foreach(irep = 1:nrep) %dorng% sample_unbiasedestimator(pcnmh_kernel, coupled_pcnmh_kernel, rinit, h = function(x) x, k = k, m = m, lag = lag)

## MSE work product = variance work product
var(sapply(rwmh_uestimators, function(l) l$uestimator[1])) * mean(sapply(rwmh_uestimators, function(l) l$cost))
var(sapply(pcnmh_uestimators, function(l) l$uestimator[1])) * mean(sapply(pcnmh_uestimators, function(l) l$cost))
## compared to 
coda::spectrum0(rwmh_chain[1e3:nrow(rwmh_chain),1])$spec
coda::spectrum0(pcnmh_chain[1e3:nrow(pcnmh_chain),1])$spec
## we lose little in terms of inefficiency if lag, k and m are chosen appropriately

## check agreement of posterior mean estimates
rowMeans(sapply(rwmh_uestimators, function(l) l$uestimator))
rowMeans(sapply(pcnmh_uestimators, function(l) l$uestimator))
colMeans(pcnmh_chain[1e3:nrow(pcnmh_chain),])

## variance of the unbiased estimator
var(sapply(rwmh_uestimators, function(l) l$uestimator[1]))
## variance of the bias correction term 
var(sapply(rwmh_uestimators, function(l) l$correction[1]))
## variance of the MCMC estimator with k-1 steps as burn-in and m steps in total 
var(sapply(rwmh_uestimators, function(l) l$mcmcestimator[1]))

