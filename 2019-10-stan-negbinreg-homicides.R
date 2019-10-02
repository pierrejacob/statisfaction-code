## This script shows how to approximate the BayesBag distribution
## introduced by Peter Buhlmann in "Discussion of Big Bayes Stories and BayesBag"
## https://www.jstor.org/stable/43288454 
## which is essentially the combination of bagging and Bayesian inference

## The example is taken from 
## https://jrnold.github.io/bugs-examples-in-stan/negbin
## and is essentially a GLM, with a Negative Binomial distribution for the outcome
## For convenience the data are included in plain text below.

## First, the data are models are defined. 
## Then plain Stan is used to approximate both the standard posterior
## and BayesBag. An issue lies in the difficulty in assessing the burn-in bias.
## Then unbiased MCMC functions are defined to provide unbiased estimators
## of expectations with respect to both the standard posterior and BayesBag.

library("tidyverse")
library("rstan")
library(unbiasedmcmc)
library(lubridate)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores())
set.seed(1)
## data in plain text
st_louis_census <- structure(list(i8094 = c(4, 5, 7, 12, 23, 56, 28, 58, 25, 50, 
                                            51, 38, 99, 71, 79, 32, 42, 46, 56, 45, 67, 33, 43, 60, 54, 43, 
                                            61, 71, 36, 40, 39, 39, 59, 59, 52, 16, 37, 56, 31, 38, 39, 93, 
                                            19, 17, 9, 37, 42, 31, 14, 28, 68, 35, 14, 34, 14, 12, 4, 2, 
                                            5, 26, 15, 4, 15, 5, 0, 60, 33, 0, 4, 19, 15, 8, 1, 3, 15, 2, 
                                            1, 8, 8, 25, 29, 7, 0, 6, 5, 10, 16, 4, 5, 18, 1, 3, 1, 15, 8, 
                                            2, 6, 10, 12, 2, 1, 2, 2, 2, 2, 2, 15, 1, 8, 14, 61), 
                                  pcunemp9 = c(16.5,8.43000030517578, 6.40000009536743, 0, 12.1999998092651, 15.4799995422363, 
                                               18.6599998474121, 15.7299995422363, 11.9899997711182, 15.5900001525879, 
                                               16.0900001525879, 17.5799999237061, 25.75, 16.0300006866455, 
                                               24.1800003051758, 10.5900001525879, 19.6399993896484, 16.4899997711182, 
                                               12.9300003051758, 18.6800003051758, 22.5, 11.1400003433228, 10.9700002670288, 
                                               17.1000003814697, 17.2299995422363, 13.1300001144409, 21.0400009155273, 
                                               24.0499992370605, 35.7599983215332, 16.9500007629395, 24.2199993133545, 
                                               20.6100006103516, 15.3100004196167, 18.9400005340576, 12.1300001144409, 
                                               10.2600002288818, 26.1700000762939, 18.8400001525879, 12, 22.1299991607666, 
                                               25.6700000762939, 15.4300003051758, 9.15999984741211, 6.5, 2.4300000667572, 
                                               10.9799995422363, 35.439998626709, 34.8499984741211, 5.15000009536743, 
                                               6.51000022888184, 24.2600002288818, 12.8000001907349, 11.960000038147, 
                                               15.960000038147, 1.22000002861023, 8.86999988555908, 11.25, 5.84999990463257, 
                                               2.79999995231628, 24.2399997711182, 11.5299997329712, 6.3899998664856, 
                                               10.8400001525879, 6.17000007629394, 6.23999977111816, 12.1899995803833, 
                                               14.1599998474121, 6.71999979019165, 6.11999988555908, 13.7399997711182, 
                                               14.6899995803833, 8.84000015258789, 2.25999999046326, 3.47000002861023, 
                                               11.6300001144409, 3.36999988555908, 4.98999977111816, 4.61999988555908, 
                                               8.19999980926514, 18.7999992370605, 14.0100002288818, 5.59999990463257, 
                                               2.51999998092651, 12.2399997711182, 5.23000001907349, 7.96000003814697, 
                                               10.1300001144409, 3.65000009536743, 4.3899998664856, 13.5600004196167, 
                                               4.53999996185303, 5.78999996185303, 2.5, 7.53999996185303, 9.55000019073486, 
                                               1.75, 5.01999998092651, 8.14000034332275, 6.59999990463257, 6.21999979019165, 
                                               5.17999982833862, 4.71999979019165, 1.91999995708466, 3.40000009536743, 
                                               7.38000011444092, 5.07999992370606, 12.289999961853, 6.51999998092651, 
                                               5.42999982833862, 7.92999982833862, 28.1599998474121), 
                                  incrs = c(16.5179996490478, 24.4790000915527, 26.121000289917, 8.42099952697754, 27.2399997711182, 
                                            25.6700000762939, 20.7609996795654, 18.1760005950928, 28.6280002593994, 
                                            19.8519992828369, 20.568000793457, 15.1680002212524, 11.2119998931885, 
                                            14.03600025177, 12.4119997024536, 19.2590007781982, 16.2900009155273, 
                                            14.5950002670288, 13.1940002441406, 14.4499998092651, 14.6750001907349, 
                                            11.4849996566772, 13.923999786377, 13.1590003967285, 14.9540004730225, 
                                            12.0570001602173, 15.121000289917, 11.8950004577637, 8.5620002746582, 
                                            16.8519992828369, 8.7480001449585, 11.6350002288818, 12.1199998855591, 
                                            9.67700004577637, 15.5120000839233, 29.6089992523193, 13.1969995498657, 
                                            10.5959997177124, 8.53600025177002, 13.5100002288818, 9.73400020599365, 
                                            12.0930004119873, 23.1650009155273, 23.7759990692139, 25.4120006561279, 
                                            14.0459995269775, 5.40700006484985, 5.44399976730347, 17.3099994659424, 
                                            11.210000038147, 7.43100023269653, 11.3839998245239, 4.99900007247925, 
                                            17.2360000610352, 24.3360004425049, 18.2019996643066, 9.3439998626709, 
                                            25.193000793457, 23.3239994049072, 9.82600021362305, 10.4169998168945, 
                                            25.6159992218018, 21.3950004577637, 20.6539993286133, 23.625, 19.5289993286133, 
                                            17.7900009155273, 22.0310001373291, 19.261999130249, 
                                            15.0340003967285, 16.7169990539551, 23.4039993286133, 26.4729995727539, 
                                            24.3099994659424, 19.5109996795654, 26.863000869751, 24.5249996185303, 
                                            23.4860000610352, 19.2280006408691, 13.53600025177, 19.9850006103516, 
                                            21.9519996643066, 31.9179992675781, 13.496000289917, 25.7700004577637, 
                                            19.375, 18.0149993896484, 20.5799999237061, 26.8479995727539, 
                                            14.9829998016357, 27.2870006561279, 21.57200050354, 32.1349983215332, 
                                            25.5839996337891, 21.128999710083, 33.548999786377, 23.6639995574951, 
                                            18.0919990539551, 16.875, 24.25, 28.1480007171631, 27.306999206543, 
                                            24.2049999237061, 33.1780014038086, 23.0620002746582, 25.6429996490478, 
                                            19.6569995880127, 28.9629993438721, 17.318000793457, 25.375, 
                                            6.27500009536743)), class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -111L))

## otherwise data are seemingly accessible via data("st_louis_census", package = "bayesjackman")

head(st_louis_census)

## create bootstrapped versions of the data
bootstrapped_datasets <- list()
nboot <- 1000
for (iboot in 1:nboot){
   bootstrapped_datasets[[iboot]] <- sample_n(st_louis_census, size = nrow(st_louis_census), replace = TRUE)
}

## function to transform data set into a list to be given to stan
transform_for_stan <- function(dataset){
   negbin_data <- within(list(), {
      y <- dataset$i8094
      N <- length(y)
      X <- model.matrix(~ 0 + pcunemp9 + incrs, data = dataset) %>% scale()
      K <- ncol(X)
      beta_mean <- rep(0, K)
      beta_scale <- rep(2.5, K)  
      alpha_mean <- 0
      alpha_scale <- 10
      reciprocal_phi_scale <- 5
   })
   return(negbin_data)
}

## Stan model of the GLM
fileName <- "negbin.stan"
write(" 
data {
  int N;
  int y[N];
  int K;
  matrix[N, K] X;
  // priors
  real alpha_mean;
  real alpha_scale;
  vector[K] beta_mean;
  vector[K] beta_scale;
  real reciprocal_phi_scale;
}
parameters {
  real alpha;
  vector[K] beta;
  real reciprocal_phi;
}
transformed parameters {
  vector[N] eta;
  real phi;
  eta = alpha + X * beta;
  phi = 1. / reciprocal_phi;
}
model {
  reciprocal_phi ~ cauchy(0., reciprocal_phi_scale);
  alpha ~ normal(alpha_mean, alpha_scale);
  beta ~ normal(beta_mean, beta_scale);
  y ~ neg_binomial_2_log(eta, phi);
}
generated quantities {
  vector[N] mu;
  vector[N] log_lik;
  vector[N] y_rep;
  mu = exp(eta);
  for (i in 1:N) {
    log_lik[i] = neg_binomial_2_log_lpmf(y[i] | eta[i], phi);
    y_rep[i] = neg_binomial_2_rng(mu[i], phi);
  }
}
", fileName)
## Compile stan model
negbin_mod <- stan_model(fileName)

## run Stan on original data set 
negbin_data <- transform_for_stan(st_louis_census)
negbin_fit <- sampling(negbin_mod, data = negbin_data, warmup = 1000, iter = 2000, chains = 4, cores = 4)
postsamples <- extract(negbin_fit)
upostsamples <- cbind(postsamples$alpha, postsamples$beta, postsamples$reciprocal_phi)

## now BayesBag, by running stan for each data set independently,
## crossing fingers that the warmup is enough in each case
## and aggregating the results

bbagsamples <- foreach (iboot = 1:nboot, .combine = rbind) %dorng% {
   standata <- transform_for_stan(bootstrapped_datasets[[iboot]])
   negbin_fit <- sampling(negbin_mod, data = standata, warmup = 1000, iter = 2000, chains = 1, cores = 1)
   postsamples <- extract(negbin_fit)
   upostsamples <- cbind(postsamples$alpha, postsamples$beta, postsamples$reciprocal_phi)
   upostsamples
}

# ## Next, compare standard posterior with BayesBag posterior
par(mfrow = c(2,2))
hist(upostsamples[,1], prob = TRUE, nclass = 100, xlim = c(2.5,3.5), main = "", xlab = expression(alpha))
hist(bbagsamples[,1], prob = TRUE, nclass = 100, add = TRUE, col = rgb(1,0,0,0.5))

hist(upostsamples[,2], prob = TRUE, nclass = 100, xlim = c(0.1,1.5), main = "", xlab = expression(beta[1]))
hist(bbagsamples[,2], prob = TRUE, nclass = 100, add = TRUE, col = rgb(1,0,0,0.5))

hist(upostsamples[,3], prob = TRUE, nclass = 100, xlim = c(-1,0.5), main = "", xlab = expression(beta[2]))
hist(bbagsamples[,3], prob = TRUE, nclass = 100, add = TRUE, col = rgb(1,0,0,0.5))

hist(upostsamples[,4], prob = TRUE, nclass = 100, xlim = c(0,1), main = "", xlab = expression(1/phi))
hist(bbagsamples[,4], prob = TRUE, nclass = 100, add = TRUE, col = rgb(1,0,0,0.5))

## The BayesBag distribution is more diffuse but otherwise overlaps quite a bit with
## the standard posterior...

## next, unbiased MCMC for both the standard posterior and BayesBag
## we will use HMC and thus we first define a bunch of tuning parameters
## HMC stepsize
stepsize <- 0.5
## HMC number of leap frog steps
nleapfrogsteps <- 5
## mass matrix
postcov <- cov(upostsamples)
postmean <- colMeans(upostsamples)
invmassdiag <- diag(cov(upostsamples))
massdiag <- 1/invmassdiag
sqrt_massdiag <- sqrt(massdiag)
## target dimension
target_dim <- get_num_upars(negbin_fit)
## function to evaluate log density on unconstrained space
stan_logtarget <- function(x) log_prob(negbin_fit, x)
## function to evaluate gradient of log density on unconstrained space
stan_gradlogtarget <- function(x) grad_log_prob(negbin_fit, x)
## stan_logtarget(upostsamples[1,]) ## test

## initial distribution of the chains
library(mvtnorm)
## this is a normal but constrained on the "logtarget" evaluated at the draw being well-defined
rinit <- function(){
   finish <- FALSE
   while (!finish){
      chain_state <- mvtnorm::rmvnorm(1, postmean, postcov)
      current_pdf <- try(stan_logtarget(chain_state))
      if (inherits(current_pdf, "try-error") || is.infinite(current_pdf) || is.na(current_pdf)){
         finish <- FALSE
      } else {
         return(list(chain_state = chain_state[1,], current_pdf = current_pdf))
      }
   }
}

## HMC kernel with fixed stepsize and fixed number of leap frog steps
hmc_kernel <- function(state){
   # draw momentum
   initial_momentum <- rnorm(target_dim, 0, sqrt_massdiag)
   chain_state <- state$chain_state
   position <- chain_state
   # leap frog integrator
   momentum <- initial_momentum + stepsize * stan_gradlogtarget(position) / 2
   for (step in 1:nleapfrogsteps){
      position <- position + stepsize * invmassdiag * momentum
      # position <- position + stepsize * (invmassmatrix %*% momentum)[,1]
      if (step != nleapfrogsteps){
         momentum <- momentum + stepsize * stan_gradlogtarget(position)
      }
   }
   momentum <- momentum + stepsize * stan_gradlogtarget(position) / 2
   proposed_pdf <- stan_logtarget(position)
   if (is.na(proposed_pdf)) proposed_pdf <- -Inf
   current_pdf <- state$current_pdf
   accept_ratio <- proposed_pdf - current_pdf
   # the acceptance ratio also features the "kinetic energy" term of the extended target
   accept_ratio <- accept_ratio + (-0.5 * sum(momentum * invmassdiag * momentum)) -
      (-0.5 * sum(initial_momentum * invmassdiag * initial_momentum))
   accept <- FALSE
   if (is.finite(accept_ratio)){
      accept <- (log(runif(1)) < accept_ratio)
   }
   if (accept){
      chain_state <- position
      current_pdf <- proposed_pdf
   }
   return(list(chain_state = chain_state, current_pdf = current_pdf, accept = accept))
}

## coupled HMC kernel, with common velocities
coupled_hmc_kernel <- function(state1, state2){
   chain_state1 <- state1$chain_state; current_pdf1 <- state1$current_pdf
   chain_state2 <- state2$chain_state; current_pdf2 <- state2$current_pdf
   # draw same momentum for two chains
   initial_momentum <- sqrt_massdiag * rnorm(target_dim, 0, 1)
   position1 <- chain_state1
   position2 <- chain_state2
   # leap frog integrator
   momentum1 <- initial_momentum + stepsize * stan_gradlogtarget(position1) / 2
   momentum2 <- initial_momentum + stepsize * stan_gradlogtarget(position2) / 2
   for (step in 1:nleapfrogsteps){
      position1 <- position1 + stepsize * invmassdiag * momentum1
      position2 <- position2 + stepsize * invmassdiag * momentum2
      if (step != nleapfrogsteps){
         momentum1 <- momentum1 + stepsize * stan_gradlogtarget(position1)
         momentum2 <- momentum2 + stepsize * stan_gradlogtarget(position2)
      }
   }
   momentum1 <- momentum1 + stepsize * stan_gradlogtarget(position1) / 2
   momentum2 <- momentum2 + stepsize * stan_gradlogtarget(position2) / 2
   proposed_pdf1 <- stan_logtarget(position1)
   proposed_pdf2 <- stan_logtarget(position2)
   if (is.na(proposed_pdf1)) proposed_pdf1 <- -Inf
   if (is.na(proposed_pdf2)) proposed_pdf2 <- -Inf
   accept_ratio1 <- proposed_pdf1 - current_pdf1
   # the acceptance ratio also features the "kinetic energy" term of the extended target
   accept_ratio1 <- accept_ratio1 + (-0.5 * sum(momentum1 * invmassdiag * momentum1)) -
      (-0.5 * sum(initial_momentum * invmassdiag * initial_momentum))
   accept_ratio2 <- proposed_pdf2 - current_pdf2
   accept_ratio2 <- accept_ratio2 + (-0.5 * sum(momentum2 * invmassdiag * momentum2)) -
      (-0.5 * sum(initial_momentum * invmassdiag * initial_momentum))
   # same uniform to accept/reject proposed state
   logu <- log(runif(1))
   accept1 <- FALSE; accept2 <- FALSE
   if (is.finite(accept_ratio1)){
      accept1 <- (logu < accept_ratio1)
   }
   if (is.finite(accept_ratio2)){
      accept2 <- (logu < accept_ratio2)
   }
   if (accept1){
      chain_state1 <- position1
      current_pdf1 <- proposed_pdf1
   }
   if (accept2){
      chain_state2 <- position2
      current_pdf2 <- proposed_pdf2
   }
   return(list(state1 = list(chain_state = chain_state1, current_pdf = current_pdf1),
               state2 = list(chain_state = chain_state2, current_pdf = current_pdf2),
               identical = FALSE))
}

## reflection maximal coupling of normal(0,D)
reflmaxcoupling <- function(mu1, mu2, sqrtD, kappa = 1){
   dim_ <- length(mu1)
   momentum1 <- rnorm(dim_, 0, 1)
   momentum2 <- rep(0, dim_)
   logu <- log(runif(1))
   z <- (mu1 - mu2) / sqrtD
   normz <- sqrt(sum(z^2))
   evector <- z / normz
   edotxi <- sum(evector * momentum1)
   if (logu < (dnorm(edotxi + kappa * normz, 0, 1, log = TRUE) - dnorm(edotxi, 0, 1, log = TRUE))){
      momentum2 <- momentum1 + kappa * z
      samesame <- TRUE
   } else {
      momentum2 <- momentum1 - 2 * edotxi * evector
      samesame <- FALSE
   }
   momentum1 <- momentum1 * sqrtD
   momentum2 <- momentum2 * sqrtD
   return(list(momentum1 = momentum1, momentum2 = momentum2, samesame = samesame))
}

## standard deviation of proposal of random walk 
sigma_mhproposal <- sqrt(invmassdiag) / 1e2
## MH kernel with random walk proposals
mh_kernel <- function(state){
   chain_state <- state$chain_state
   current_pdf <- state$current_pdf
   proposal_value <- chain_state + sigma_mhproposal * rnorm(target_dim)
   proposal_pdf <- stan_logtarget(proposal_value)
   accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
   if (accept){
      return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
   } else {
      return(list(chain_state = chain_state, current_pdf = current_pdf))
   }
}
## coupled MH kernel, with "reflection-maximally" coupled proposals
coupled_mh_kernel <- function(state1, state2){
   chain_state1 <- state1$chain_state; current_pdf1 <- state1$current_pdf
   chain_state2 <- state2$chain_state; current_pdf2 <- state2$current_pdf
   proposal_value <- reflmaxcoupling(chain_state1, chain_state2, sigma_mhproposal, kappa = 1)
   proposal1 <- chain_state1 + proposal_value$momentum1
   proposal_pdf1 <- stan_logtarget(proposal1)
   if (proposal_value$samesame){
      proposal2 <- proposal1; proposal_pdf2 <- proposal_pdf1
   } else {
      proposal2 <- chain_state2 + proposal_value$momentum2
      proposal_pdf2 <- stan_logtarget(proposal2)
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
   identical_ <- proposal_value$samesame && accept1 && accept2
   return(list(state1 = list(chain_state = chain_state1, current_pdf = current_pdf1),
               state2 = list(chain_state = chain_state2, current_pdf = current_pdf2),
               identical = identical_))
}
## mixture kernels, equal to MH with probability omega and HMC otherwise
omega <- 1/10
mixture_single_kernel <- function(state){
   if (runif(1) < omega){
      return(mh_kernel(state))
   } else {
      return(hmc_kernel(state))
   }
}
## coupled mixture kernel
mixture_coupled_kernel <- function(state1, state2){
   if (runif(1) < omega){
      return(coupled_mh_kernel(state1, state2))
   } else {
      return(coupled_hmc_kernel(state1, state2))
   }
}

## lag between the chains
lag <- 1e3
## sample meeting times
# sample_meetingtime(mixture_single_kernel, mixture_coupled_kernel, rinit, lag = lag)
nrep <- 1e2
meetingtimes <- foreach(irep = 1:nrep, .combine = c) %dorng% {
   sample_meetingtime(mixture_single_kernel, mixture_coupled_kernel, rinit, lag = lag)$meetingtime
}
## Histogram of times to meet, counting from the time the chains start being coupled
par(mfrow = c(2,1))
hist(meetingtimes - lag, prob = TRUE)
# TV upper bounds as in Biswas, Jacob & Vanetti 2019 https://arxiv.org/abs/1905.09971
tv_upper_bound_estimates <- function(coupling_times, L, t){
   return(pmax(0,ceiling((coupling_times-L-t)/L)))
}
## plot TV upper bounds as a function of iteration
timerange <- floor(seq(from = 1, to = 1000, length.out = 100))
tvupperbounds <- sapply(timerange, function(t) mean(tv_upper_bound_estimates(meetingtimes, lag, t)))
plot(timerange, tvupperbounds, type = "l",
     xlab = "iteration", ylab = "TV upper bounds")
## this indicates that warm-up takes a few hundred steps in this case

## unbiased MCMC approximation of standard posterior 
nrep <- 1000
# length of pairs of chains
m <- 2e3 #
coupled_chains <- foreach(irep = 1:nrep) %dorng% {
   sample_coupled_chains(mixture_single_kernel, mixture_coupled_kernel, rinit, m = m, lag = lag)
}
## now let's see if unbiased MCMC agrees with stan on the standard posterior
hist1 <- histogram_c_chains(coupled_chains, 1, lag, m)
hist1 <- plot_histogram(hist1) + geom_density(data = data.frame(x = upostsamples[,1]), aes(x = x, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL))
hist1 <- hist1 + xlab(expression(alpha)) + xlim(2.5, 3.4)
#
hist2 <- histogram_c_chains(coupled_chains, 2, lag, m)
hist2 <- plot_histogram(hist2) + geom_density(data = data.frame(x = upostsamples[,2]), aes(x = x, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL))
hist2 <- hist2 + xlab(expression(beta[1])) + xlim(0, 1.4)
#
hist3 <- histogram_c_chains(coupled_chains, 3, lag, m)
hist3 <- plot_histogram(hist3) + geom_density(data = data.frame(x = upostsamples[,3]), aes(x = x, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL))
hist3 <- hist3 + xlab(expression(beta[2])) + xlim(-1, 0.25)
#
hist4 <- histogram_c_chains(coupled_chains, 4, lag, m)
hist4 <- plot_histogram(hist4) + geom_density(data = data.frame(x = upostsamples[,4]), aes(x = x, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL))
hist4 <- hist4 + xlab(expression(1/phi)) + xlim(0, 1)
gridExtra::grid.arrange(hist1, hist2, hist3, hist4, nrow = 2, ncol = 2)
## seems to agree

## unbiased MCMC approximation of the BayesBag posterior
## basically the same as above, but the data is bootstrapped before running the coupled chains
coupled_chains_bayesbag <- foreach(irep = 1:nrep) %dorng% {
   standata <- transform_for_stan(bootstrapped_datasets[[irep]])
   # very quick call to "sampling" to have access to the target functions
   # ... (seems a bit hacky)...
   negbin_fit <- sampling(negbin_mod, data = standata, warmup = 1, iter = 2, chains = 1, cores = 1)
   ## to make sure the following definitions overwrite the previous definitions 
   ## we can use double arrows... seems a bit hacky 
   ## function to evaluate log density on unconstrained space
   stan_logtarget <<- function(x) log_prob(negbin_fit, x)
   ## function to evaluate gradient of log density on unconstrained space
   stan_gradlogtarget <<- function(x) grad_log_prob(negbin_fit, x)
   ## run coupled chains
   sample_coupled_chains(mixture_single_kernel, mixture_coupled_kernel, rinit, m = m, lag = lag)
}
summary(sapply(coupled_chains_bayesbag, function(x) x$meetingtime))
summary(sapply(coupled_chains_bayesbag, function(x) x$cost))

## see agreement between unbiased MCMC and stan on BayesBag
k <- lag
hist1bb <- histogram_c_chains(coupled_chains_bayesbag, 1, k, m)
hist1bb <- plot_histogram(hist1bb) + geom_density(data = data.frame(x = bbagsamples[,1]), aes(x = x, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL))
hist1bb <- hist1bb + xlab(expression(alpha)) + xlim(2.5, 3.4)
# 
hist2bb <- histogram_c_chains(coupled_chains_bayesbag, 2, k, m)
hist2bb <- plot_histogram(hist2bb) + geom_density(data = data.frame(x = bbagsamples[,2]), aes(x = x, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL))
hist2bb <- hist2bb + xlab(expression(beta[1])) + xlim(0, 1.4)
#
hist3bb <- histogram_c_chains(coupled_chains_bayesbag, 3, k, m)
hist3bb <- plot_histogram(hist3bb) + geom_density(data = data.frame(x = bbagsamples[,3]), aes(x = x, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL))
hist3bb <- hist3bb + xlab(expression(beta[2])) + xlim(-1, 0.25)
#
hist4bb <- histogram_c_chains(coupled_chains_bayesbag, 4, k, m)
hist4bb <- plot_histogram(hist4bb) + geom_density(data = data.frame(x = bbagsamples[,4]), aes(x = x, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL))
hist4bb <- hist4bb + xlab(expression(1/phi)) + xlim(0, 1)

gridExtra::grid.arrange(hist1bb, hist2bb, hist3bb, hist4bb, nrow = 2, ncol = 2)
## seems to agree!

## compare cumulative distribution functions (CDFs) of standard Bayes and BayesBag
## for parameters beta[1] and beta[2]

## obtain signed measure representations of the distributions
cchainsdf_ <- unbiasedmcmc::c_chains_to_dataframe(coupled_chains, k, m, dopar = TRUE)
cchainsbbdf_ <- unbiasedmcmc::c_chains_to_dataframe(coupled_chains_bayesbag, k, m, dopar = TRUE)
subiterations <- floor(seq(from = 1, to = nrow(cchainsdf_2), length.out = 1000))
subiterationsbb <- floor(seq(from = 1, to = nrow(cchainsbbdf_2), length.out = 1000))
## sort by ascending order of certain components and sum weights
cchainsdf_2 <- cchainsdf_ %>% select(rep, weight, atom.2) %>% arrange(atom.2) %>% mutate(cumw = cumsum(weight))
cchainsbbdf_2 <- cchainsbbdf_ %>% select(rep, weight, atom.2) %>% arrange(atom.2) %>% mutate(cumw = cumsum(weight))
cchainsdf_3 <- cchainsdf_ %>% select(rep, weight, atom.3) %>% arrange(atom.3) %>% mutate(cumw = cumsum(weight))
cchainsbbdf_3 <- cchainsbbdf_ %>% select(rep, weight, atom.3) %>% arrange(atom.3) %>% mutate(cumw = cumsum(weight))

## plot approximations of CDFs
par(mfrow = c(1,2))
plot(x = cchainsdf_2[subiterations,]$atom.2, y = cchainsdf_2[subiterations,]$cumw, type = "l", col = "blue", 
     xlab = expression(beta[1]), ylab = "CDF")
lines(x = cchainsbbdf_2[subiterationsbb,]$atom.2, y = cchainsbbdf_2[subiterationsbb,]$cumw, lwd = 2, col = "red")

plot(x = cchainsdf_3[subiterations,]$atom.3, y = cchainsdf_3[subiterations,]$cumw, type = "l", col = "blue",
     xlab = expression(beta[2]), ylab = "CDF")
lines(x = cchainsbbdf_3[subiterationsbb,]$atom.3, y = cchainsbbdf_3[subiterationsbb,]$cumw, lwd = 2, col = "red")


