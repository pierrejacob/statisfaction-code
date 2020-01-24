library(reshape2)
library(ggplot2)
theme_set(theme_classic())
theme_update(axis.text.x = element_text(size = 20),
             axis.text.y = element_text(size = 20),
             axis.title.x = element_text(size = 25, margin=margin(20,0,0,0)),
             axis.title.y = element_text(size = 25, angle = 90, margin = margin(0,20,0,0)),
             strip.text = element_text(size = 25),
             strip.background = element_rect(fill="white"),
             legend.position = "bottom")

set.seed(1)
rm(list = ls())

## target distribution of the MCMC algorithm
library(Rcpp)
cppFunction("
double target_cpp(double x){
  return log(abs(cos(x)) * (exp(-0.5 * (x + 2.) * (x + 2.)) + exp(-0.5 * (x - 2.) * (x - 2.))));
}")
target <- target_cpp

## initial distribution of the chain
initmean <- 20
initsd <- 1
rinit <- function(){
  chain_state <- rnorm(1, initmean, initsd)
  current_pdf <- target(chain_state)
  return(list(chain_state = chain_state, current_pdf = current_pdf))
}

## Metropolis-Hastings transition kernel
sd_proposal <- 1.5
single_kernel <- function(state){
  proposal_value <- rnorm(1, mean=state$chain_state, sd=sd_proposal)
  proposal_pdf <- target(proposal_value)
  accept <- (log(runif(1)) < (proposal_pdf - state$current_pdf))
  if (accept){
    return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
  } else {
    return(list(chain_state = state$chain_state, current_pdf = state$current_pdf))
  }
}

## number of MCMC iterations
nmcmc <- 300
## generate Markov chain
current <- rinit()
chain <- rep(0, nmcmc)
for (imcmc in 1:nmcmc){
  current <- single_kernel(current)
  chain[imcmc] <- current$chain_state
}

traceplot <- qplot(x = 1:nmcmc, y = chain, geom = "line") + xlab('iteration') + ylab('x')
traceplot

## Now pretend that we forgot the value of 'sd_proposal'
## Our goal is to estimate it based on the MCMC trace (in 'chain')

## Denote by alpha(x,x') the MH acceptance ratio, here equal to min(1, pi(x')/pi(x))
## and denote by a(x) the probability of accepting any proposal from x,
## a(x) = integral of q(x,x') alpha(x,x') dx'

## The MH transition kernel has a density
## K(x, x') =
##            q(x,x') alpha(x,x') if x is different from x'
##            1 - a(x)            if x is equal to x'
## with respect to Lebesgue + Dirac on x

## The hard part is that a(x) is not analytically available.

## We consider the following estimator of 1-a(x), i.e. of the rejection probability from x
estim_rejectionproba_basic <- function(x, sigma, size = 1){
  proposal_value <- rnorm(size, mean=x, sd=sigma)
  mh_ratios <- sapply(proposal_value, function(v) pmin(1, exp(target(v) - target(x))))
  return(mean(1-mh_ratios))
}
## The problem with the above estimator is that it can be equal to zero exactly.

## We consider the following alternative estimator of 1-a(x).
## The justification is based on the fact that if X is Negative Binomial with parameters r,p
## then Neuts and Zacks (1967) state that E[(r + X − 1)−1] = (1 − p)/(r − 1).
## This is explained e.g. in "The Alive Particle Filter" by Jasra, Lee, Yau, Zhang https://arxiv.org/abs/1304.0151
estim_rejectionproba_sophisticated <- function(x, sigma, size = 2){
  success <- 0
  ntrials <- 0
  while (success < size){
    ntrials <- ntrials + 1
    proposal_value <- rnorm(1, mean=x, sd=sigma)
    mh_ratio <- min(1, exp(target(proposal_value) - target(x)))
    success <- success + (runif(1) < (1-mh_ratio))
  }
  return((size - 1) / (ntrials-1))
}

## We can also compute the rejection probability by numerical integration
## in this very simple example.
# NI_transition <- function(x, sigma){
#   targetx <- target(x)
#   mh_ratio <- function(z) sapply(z, function(v) min(1, exp(target(v) - targetx)))
#   f_ <- function(z) mh_ratio(z) * dnorm(z, mean = x, sd = sigma)
#   return(1 - integrate(f_, lower = -Inf, upper = +Inf, subdivisions = 1e4)$value)
# }

## Now we can test the above functions
# idxsame <- which(diff(chain)==0)[1]
# chain[idxsame] == chain[idxsame+1]
# x <- chain[idxsame]
# sigma <- 3
# NI_transition(x, sigma)
# estim_rejectionproba_basic(x, sigma = sigma, size = 1e5)
# estim_rejectionproba_sophisticated(x, sigma = sigma, size = 1e5)

## likelihood estimator
get_ll <- function(sigma_seq, size, type = "basic"){
  ll <- rep(0, length(sigma_seq))
  for (i in 1:length(sigma_seq)){
    sigma <- sigma_seq[i]
    ll[i] <- 0
    ## go through chain
    for (it in 2:nmcmc){
      ## either state is same as previous one or it's different
      if (identical(chain[it], chain[it-1])){
        ## if identical, sample likelihood estimator
        if (type == "basic"){
          ll[i] <- ll[i] + log(estim_rejectionproba_basic(chain[it-1], sigma = sigma, size = size))
        } else {
          ll[i] <- ll[i] + log(estim_rejectionproba_sophisticated(chain[it-1], sigma = sigma, size = size))
        }
      } else {
        ll[i] <- ll[i] + dnorm(chain[it], chain[it-1], sd = sigma, log = TRUE) + log(min(1, exp(target(chain[it]) - target(chain[it-1]))))
      }
    }
    ll[i]
  }
  return(ll)
}

## Test of the log-likelihood estimator
size <- 2
sigma_seq <- seq(from = 1, to = 5, length.out = 20)
get_ll(sigma_seq[1], size = 2, type = "basic")
get_ll(sigma_seq[1], size = 2, type = "nonbasic")
## The "basic" estimator will often by zero whereas the non-basic one will never be zero.

nrep <- 100
size <- 2
##
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)

llestimators_nonbasic <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  get_ll(sigma_seq, size, type = "nonbasic")
}

llestimators_nonbasic_df <- reshape2::melt(llestimators_nonbasic)
names(llestimators_nonbasic_df) <- c("rep", "isigma", "value")
sigma_df <- data.frame(isigma = 1:length(sigma_seq), sigma = sigma_seq)
llestimators_nonbasic_df <- merge(llestimators_nonbasic_df, sigma_df, by =  "isigma")

ggplot(llestimators_nonbasic_df, aes(x = sigma, y = value, group = isigma)) + geom_boxplot() + xlab(expression(sigma)) + ylab("log-likelihood")

### Now we can envision Bayesian inference on sigma
### prior: Uniform on [1,5]
logdprior <- function(x) dunif(x, min = 1, max = 5, log = TRUE)
## pseudo-marginal algorithm
niterations <- 1e3
size <- 50
sigma_chain <- rep(2, niterations)
target_current <- logdprior(sigma_chain[1]) + get_ll(sigma_seq = sigma_chain[1], size = size, type = "nonbasic")
for (iter in 2:niterations){
  print(iter)
  sigma_proposal <- rnorm(n = 1, mean = sigma_chain[iter-1], sd = .1)
  target_proposal <- logdprior(sigma_proposal)
  if (is.finite(target_proposal)){
    target_proposal <- target_proposal + get_ll(sigma_seq = sigma_proposal, size = size, type = "nonbasic")
  }
  u <- runif(1)
  if (log(u) < (target_proposal - target_current)){
    target_current <- target_proposal
    sigma_chain[iter] <- sigma_proposal
  } else {
    sigma_chain[iter] <- sigma_chain[iter-1]
  }
}

##
gmeta <- qplot(x = 1:length(sigma_chain), y = sigma_chain, geom = "blank") + geom_line() + xlab("iteration") + ylab(expression(sigma))
gmeta
