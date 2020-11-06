# Custom survival distributions in R
library(survival)
library(flexsurv)
library(ggplot2)
library(cowplot)
library(dplyr)
library(magrittr)

# https://stats.stackexchange.com/questions/357621/how-is-the-weibull-distribution-parameterised-in-survreg-in-r

n <- 150
nsim <- 250

################################################
# ---- EXAMPLE 1: GETTING WEIBULL TO WORK ---- #

# rweibullPH parameterizes the shape (a) and scale (m) as follows
# f(t) = a * m * t^(a-1) * exp(-m*t^a)
# When shape=a=1, then we have the exponential

shape <- 0.75
scale <- 2.25
holder <- matrix(nrow=nsim,ncol=2)
set.seed(1234)
for (ii in seq(nsim)) {
  t2e <- flexsurv::rweibullPH(n=n, shape=shape, scale=scale)
  dd <- rep(1,n)
  df <- data.frame(t2e=t2e,dd=dd)
  fit <- flexsurvreg(Surv(t2e,dd)~1,data=df,dist='weibull')$coefficients
  shape_trans <- exp(fit['shape'])
  scale_trans <- (1/exp(fit['scale']))^shape_trans
  holder[ii,] <- c(shape_trans, scale_trans)
}
apply(holder,2,mean)

# Function to recover quantile
qweibullPH(p=0.5, shape=shape_trans, scale=scale_trans)
1-exp(-scale_trans*(median(t2e)^shape_trans))
median(t2e)
(log(2)/scale_trans)^(1/shape_trans)
(1/scale_trans^(1/shape_trans))*log(2)^(1/shape_trans)

mu_weibullPH <- function(shape, scale) {
  return( gamma(1+1/shape) / (scale^(1/shape)) )
}

mu_weibullPH(shape=shape_trans, scale=scale_trans)
mean(t2e)

#############################################
# ---- EXAMPLE 2: Multivariate Weibull ---- #

shape <- 0.75
b0 <- c(1, -1, 1)
p <- length(b0)
holder <- matrix(nrow=nsim,ncol=p+1)
set.seed(1234)
for (ii in seq(nsim)) {
  X <- cbind(1,matrix(runif(n*(p-1)),ncol=p-1)) %>% set_colnames(c('int','X1','X2'))
  eta <- as.vector(X %*% b0)
  t2e <- rweibullPH(n=n, shape=shape, scale=exp(eta))
  dd <- rep(1,n)
  df <- cbind(data.frame(t2e=t2e,dd=dd),X[,-1])
  fit <- flexsurvreg(Surv(t2e,dd)~X1+X2,data=df,dist='weibull')$coefficients
  (shape_trans <- exp(fit['shape']))
  (scale_trans <- log((1/exp(fit[2]))^shape_trans))
  (bvec <- log((1/exp(fit[3:4]))^shape_trans))
  holder[ii,] <- c(shape_trans, scale_trans, bvec)
}
apply(holder, 2, mean)

# Confirm that average of averages is about the average....
mean(t2e)
sapply(exp(as.vector(X %*% c(scale_trans,bvec))), function(ss) mu_weibullPH(shape=shape_trans, scale=ss)) %>% mean


####################################################
# ---- EXAMPLE 3: GETTING EXPONENTIAL TO WORK ---- #

rate <- 0.7
holder <- matrix(nrow=nsim,ncol=1)
set.seed(1234)
for (ii in seq(nsim)) {
  t2e <- rexp(n=n, rate=rate)
  dd <- rep(1,n)
  df <- data.frame(t2e=t2e,dd=dd)
  fit <- flexsurvreg(Surv(t2e,dd)~1,data=df,dist='exponential')$coefficients
  rate_trans <- exp(fit)
  holder[ii] <- rate_trans
}
apply(holder,2,mean)

# Check median
qexp(p=0.5,rate=rate_trans)
median(t2e)

# Check mean
1/rate_trans
mean(t2e)

#################################################
# ---- EXAMPLE 4: Multivariate Exponential ---- #

b0 <- c(0.75, -1.75, 1.75)
p <- length(b0)
holder <- matrix(nrow=nsim,ncol=p)
set.seed(1234)
for (ii in seq(nsim)) {
  X <- cbind(1,matrix(runif(n*(p-1)),ncol=p-1)) %>% set_colnames(c('int','X1','X2'))
  eta <- as.vector(X %*% b0)
  t2e <- rexp(n=n, rate=exp(eta))
  dd <- rep(1,n)
  df <- cbind(data.frame(t2e=t2e,dd=dd),X[,-1])
  fit <- flexsurvreg(Surv(t2e,dd)~X1+X2,data=df,dist='exponential')$coefficients
  bvec <- fit
  holder[ii,] <- bvec
}
apply(holder, 2, mean)

# Confirm that average of averages is about the average....
mean(t2e)
sapply(exp(as.vector(X %*% bvec)), function(ss) 1/ss) %>% mean



