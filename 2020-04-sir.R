## script accompanying blog post on Statisfation in April 2020
## about Kermack and McKendrick 1927 paper

## load packages
library(deSolve)
library(ggplot2)
library(dplyr)

# we will also need reshape2, ggthemes
rm(list = ls())
## tune graphical settings
theme_set(ggthemes::theme_tufte(ticks = TRUE))
theme_update(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), 
             axis.title.x = element_text(size = 25, margin = margin(20, 0, 0, 0), hjust = 1), 
             axis.title.y = element_text(size = 25, angle = 0, margin = margin(0, 20, 0, 0), vjust = 1), 
             legend.text = element_text(size = 20), 
             legend.title = element_text(size = 20), title = element_text(size = 30), 
             strip.text = element_text(size = 25), strip.background = element_rect(fill = "white"), 
             legend.position = "bottom")

## Data for example of the article, Bombay plague 1905 1906 
## e.g. https://github.com/bbolker/fitsir/tree/master/data

plaguedatatxt <- 
  " 1 5
  2 10
  3 17
  4 22
  5 30
  6 50
  7 51
  8 90
  9 120
  10 180
  11 292
  12 395
  13 445
  14 775
  15 780
  16 700
  17 698
  18 880
  19 925
  20 800
  21 578
  22 400
  23 350
  24 202
  25 105
  26 65
  27 55
  28 40
  29 30
  30 20"

plaguedata <- read.table(text = plaguedatatxt)
## plot data
gplague <- ggplot(data = plaguedata, aes(x = V1, y = V2)) + geom_point()
gplague <- gplague + xlab("weeks") + ylab("") + labs(subtitle = "# deaths per week", title = "plague in Bombay 1905-1906")
ggplague <- gplague + scale_x_continuous(breaks = c(5,10,15,20,25,30))
gplague <- gplague + scale_y_continuous(breaks = 1:9 * 100)
gplague <- gplague + ggthemes::geom_rangeframe()
gplague
# ggsave(filename = "plaguebombay.png", plot = gplague, width = 7, height = 7)

## We also consider another data set 
## Boarding school influenza oubteak 1976
## can be loaded e.g. with
library(epimdr)
data(flu)
## plot data
gflu <- ggplot(flu, aes(x = day, y = cases)) + geom_point()
gflu <- gflu + xlab("days") + ylab("") + labs(subtitle = "# confined in bed", title = "Influenza in a boarding school")
gflu <- gflu + scale_x_continuous(breaks = 1:14)
gflu <- gflu + scale_y_continuous(breaks = 0:3 * 100)
gflu <- gflu + ggthemes::geom_rangeframe()
gflu
# ggsave(filename = "fluboarding.png", plot = gflu, width = 7, height = 7)

## solve ODE describing SIR model using deSolve
## assumes 1 infected individual at time zero
SIR.model <- function(t, ell, kappa, N, by = 0.1){
  init_SIR <- c(S=N-1, I=1, R=0)
  parms <- c(ell=ell, kappa=kappa)
  times <- seq(from = 0, to = t, by= by)
  equations <- function(time,state,parms){
    with(as.list(c(state,parms)),{
      dS <- -kappa*S*I
      dI <- +kappa*S*I -ell*I
      dR <- +ell*I
      return(list(c(dS,dI,dR)))})
  }
  out <- ode(y=init_SIR, times=times, equations, parms=parms)
  out.df<-as.data.frame(out)
  return(out.df)
}

## show some solutions of SIR equations
ggplot(reshape2::melt(SIR.model(100, 0.15, 4e-4, 1e3, by = 1), "time"),
       aes(x=time, y=value, colour=variable)) + geom_line() + ylab("") + viridis::scale_colour_viridis(name = "", discrete=T)

## following code that is commented out
## produces the animation in the blog post 
## it's commented because it takes ~ 1 minute to run 
## and requires gganimate

# library(gganimate)
# ## create animation as N varies
# dfs <- lapply(floor(seq(from = 200, to = 1200, length.out = 15)), function(N){
#   df <- reshape2::melt(SIR.model(100, 0.15, 4e-4, N, by = 1), "time")
#   df$N <- N
#   df
# })
# dfs <- do.call(rbind, dfs)
# g <- ggplot(dfs, aes(x=time, y=value/N, colour=variable, linetype = variable)) 
# g <- g + geom_line() + ylab("") + scale_color_manual(name = "", values = c(rgb(1,0.8,0.2), 
#                                                                            rgb(.8,0.1,0.6), 
#                                                                            rgb(.15,0.5,0.8)))
# g <- g + scale_linetype(name = "")
# g <- g + transition_time(N)
# g <- g + ggtitle(label = "", subtitle = "N = {floor(frame_time)}")
# # g <- g + labs(title = "", subtitle = 'N = {closest_state}')
# # g <- g + transition_states(N, transition_length = 3, state_length = 0) + ease_aes('cubic-in-out')
# 
# animate(g, fps = 10, width = 450, height = 450)
# anim_save("anim.gif", animation = last_animation())
##

## Next we can try to reproduce the figure Kermack and McKendrick
dzdt <- function(x) 890 * (1/cosh(0.2*x-3.4))^2
xseq <- c(5.0,  7.5, 10.0, 12, 14, 15, 17, 19, 20.0, 22, 23, 25.0, 27, 30.0)
yseq <- dzdt(xseq)

gplague2 <- gplague + stat_function(fun = dzdt) + xlab("weeks") + ylab("")
gplague2 <- gplague2 + geom_point(data = data.frame(x = xseq, y = yseq), aes(x = x, y = y), size = 4, colour = 'black')
gplague2 <- gplague2 + geom_point(data = data.frame(x = xseq, y = yseq), aes(x = x, y = y), size = 3, colour = 'white')
gplague2 <- gplague2 + geom_point(data = data.frame(x = xseq, y = yseq), aes(x = x, y = y), size = .3, colour = 'black')
gplague2 <- gplague2 + ggthemes::geom_rangeframe()
gplague2
## pretty close!

## we can also try to optimize a bit 
## try to find the curve dz/dt = A / (cosh^2(Bt - phi)) with parameters (A,B,phi)
dzdt_error <- function(par){
  times <- plaguedata$V1
  fittedvalues <- sapply(times, function(x) par[1] / (cosh(par[2]*x - par[3]))^2)
  errors <- plaguedata$V2 - fittedvalues
  return(sum(errors^2))
}

optim_results <- optim(c(890, .2, 3), dzdt_error)
par <- optim_results$par
times <- plaguedata$V1
fittedvalues <- sapply(times, function(x) par[1] / (cosh(par[2]*x - par[3])^2))
new_dzdt <- function(x) par[1] * (1/cosh(par[2]*x-par[3]))^2
gplague2 + stat_function(fun = new_dzdt, colour = "red", linetype = 2) 
## new curve is very close to the previous one, a bit shifted towards the right 

## next we can try to fit the SIR model itself, instead of the above approximation
## squared error function to optimize
sir_error <- function(par){
  logell <- par[1]
  logkappa <- par[2]
  logN <- par[3]
  output <- SIR.model(30, exp(logell), exp(logkappa), exp(logN), by = 1)
  dzdt <- diff(output$R)
  return(sum((dzdt - plaguedata$V2)^2))
}

## optimize starting from with these vlaues
N <- 60000
kappa <- 1e-4
ell <- 4.32
par <- c(log(ell), log(kappa), log(N))
optim_results <- optim(par, sir_error, method = 'Nelder-Mead')
exp(optim_results$par)

SIR_output <- SIR.model(30, 5.68, 8.025e-5, 75500, by = 1)
SIR_diffR <- data.frame(time = 1:30, diffR = diff(SIR_output$R))
gplague2 + geom_line(data=SIR_diffR, aes(x = time, y = diffR), colour = 'red', linetype = 2)
## very close as well

## back to the flu data set 
N <- 763
sir_error <- function(par){
  logell <- par[1]
  logkappa <- par[2]
  output <- SIR.model(t = nrow(flu), ell = exp(logell), kappa = exp(logkappa), N = N, by = 1)
  return(sum((output$I[2:(nrow(flu)+1)] - flu$cases)^2))
}
par <- c(log(0.5), log(2/N))
optim_results <- optim(par, sir_error, method = 'Nelder-Mead')
exp(optim_results$par)

SIR_output <- SIR.model(nrow(flu), .44, 2.2e-3, N, by = .01)
gflu + geom_line(data=SIR_output, aes(x = time, y = I), colour = 'red', linetype = 2)
## we see a good fit with these parameter values

## R0 for this particular outbreak would be calculated as N * 2.2e-3 / 0.44



