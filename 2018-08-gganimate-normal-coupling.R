library(ggplot2)
library(gganimate)
set.seed(1)

# normal means
mu1 <- 0
mu2 <- 2
# std deviation
sigma <- 1
# number of samples
nsamples <- 1000
# reflection-maximal coupling
reflmax_samples <- matrix(nrow = nsamples, ncol = 2)
# draw x components first
reflmax_samples[,1] <- rnorm(nsamples)
# this follows the notation of Bou Rabee et al, 2018, roughly
z <- mu1 - mu2
normz <- sqrt(sum(z^2))
e <- z / normz
utilde <- runif(nsamples, 0, 1)
accepts <- (log(utilde) < (dnorm(e * reflmax_samples[,1] + normz, 0, 1, log = TRUE) - dnorm(e*reflmax_samples[,1], log = TRUE)))
eta <- reflmax_samples[,1]
eta[accepts] <- reflmax_samples[accepts,1] + z
eta[!accepts]<- reflmax_samples[!accepts,1] - 2 * (e * reflmax_samples[!accepts,1]) * e
reflmax_samples[,2] <- mu2 + eta

df <- data.frame(coupling = rep("reflection-maximal", nsamples), x = reflmax_samples[,1], y = reflmax_samples[,2])

# reflection coupling
refl_samples <- matrix(0, nrow = nsamples, ncol = 2)
refl_samples[,1] <- reflmax_samples[,1]
refl_samples[,2] <- (mu2-mu1)-refl_samples[,1]

df <- rbind(df, data.frame(coupling = rep("reflection", nsamples), x = refl_samples[,1], y = refl_samples[,2]))

# optimal transport coupling
transport_samples <- matrix(0, nrow = nsamples, ncol = 2)
transport_samples[,1] <- reflmax_samples[,1]
transport_samples[,2] <- mu2 - mu1 + transport_samples[,1]

df <- rbind(df, data.frame(coupling = rep("optimal transport", nsamples), x = transport_samples[,1], y = transport_samples[,2]))

# max coupling
max_samples <- matrix(0, nrow = nsamples, ncol = 2)
max_samples[,1] <- reflmax_samples[,1]
dp <- function(x) dnorm(x, mean = mu1, sd = 1, log = TRUE)
dq <- function(x) dnorm(x, mean = mu2, sd = 1, log = TRUE)
rq <- function(n) rnorm(n, mean = mu2, sd = 1)
for (isample in 1:nsamples){
  x <- max_samples[isample,1]
  if (dp(x) + log(runif(1)) < dq(x)){
    max_samples[isample,2] <- x
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rq(1)
      reject <- (dq(y) + log(runif(1)) < dp(y))
    }
    max_samples[isample,2] <- y
  }
}
df <- rbind(df, data.frame(coupling = rep("maximal", nsamples), x = max_samples[,1], y = max_samples[,2]))

## Scatter plots and marginals
ggplot(df, aes(x=x, y=y, group = coupling, colour = factor(coupling))) +
  geom_point()+
  theme_minimal() + viridis::scale_color_viridis(discrete=T)
#
ggplot(df, aes(x=x, group = coupling, fill = factor(coupling))) + geom_histogram(aes(y = ..density..), position = position_dodge()) +
  theme_minimal() + viridis::scale_fill_viridis(discrete=T) + stat_function(fun = function(x) dnorm(x, mu1, sigma))
#
ggplot(df, aes(x=y, group = coupling, fill = factor(coupling))) + geom_histogram(aes(y = ..density..), position = position_dodge()) +
  theme_minimal() + viridis::scale_fill_viridis(discrete=T) + stat_function(fun = function(x) dnorm(x, mu2, sigma))

## gganimate
ggplot(df, aes(x = x, y = y)) + xlim(-4, 4) + ylim(-5, 6) +
  geom_point()+ geom_text(data = data.frame(coupling = unique(df$coupling)), aes(label = coupling, x = -1, y = -4), size = 10) +
  theme_minimal() +
  transition_states(coupling, 3, 1)

