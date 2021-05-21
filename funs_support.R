# Utility functions

# EMULATE DEFAULT GGPLOT
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# FUNCTION WRAPPER TO GET SURV CIs from SIMULATION
# bhat: CoxPh coefficient
# Eta: exp(X bhat)
# Sigma: Var(bhat)
# Y: Surv(time,event)
# x: patient level data
pm_surv = function(bhat, Sigma, Eta, Y, x, nsim=100, alpha=0.05) {
  set.seed(nsim)
  if (!is.matrix(x)) { x = as.matrix(x) }
  stopifnot(length(Eta) == length(Y))
  stopifnot(ncol(x) == length(bhat))
  # Different beta
  bhat_sim = rmvnorm(n=nsim, mean=bhat, sigma=Sigma)
  # Needs to be the right dimensions, values are irrelevant
  X2 = matrix(rep(0,length(Y)*ncol(x)),ncol=ncol(x))
  holder = list()
  for (i in seq(nsim)) {
    bhat_i = bhat_sim[i,]
    eta_i = exp(x %*% bhat_i)[1,1]
    res_i = coxsurv.fit(ctype=1, stype = 1, se.fit = FALSE, cluster = NULL,
                        varmat = Sigma,
                        y = Y, risk = Eta, risk2 = eta_i, 
                        x = X2,
                        wt = NULL, y2 = NULL, x2 = NULL,
                        position = NULL, strata = NULL, oldid = NULL,
                        strata2 = NULL, id2 = NULL,unlist = TRUE)
    # res_i = coxsurv.fit(ctype=1, stype = 1, se.fit = FALSE, cluster = NULL,
    #                     varmat = Sigma, y = Y, x = X, wt = rep(1,nrow(X)),
    #                     risk = Eta, y2 = y, x2 = x, risk2 = eta_i,
    #                     position = NULL, strata = NULL, oldid = NULL,
    #                     strata2 = NULL, id2 = NULL,unlist = TRUE)
    res_i = mutate(as.data.frame(do.call('cbind',res_i[c('time','surv')])),idx=i)
    holder[[i]] = res_i
  }
  res_sim = as_tibble(do.call('rbind',holder))
  res_sim = res_sim %>% group_by(time) %>%  
    summarise(mu=mean(surv),lb=quantile(surv,alpha/2),ub=quantile(surv,1-alpha/2)) %>%
    arrange(time)
  return(res_sim)
}
