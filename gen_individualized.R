# --- SCRIPT TO GENERATE INDIVIDUALIZED SURVIVAL CURVES --- #

# 'magrittr','stringr','forcats', 'readr',
pckgs <- c('dplyr','tibble',
           'cowplot','ggplot2',
           'survival','mvtnorm')
for (pp in pckgs) { library(pp,character.only = T)}

source('funs_support.R')

# --- (1) LOAD MODEL --- #
load('cox_mdl.RData')

# --- (2) PATIENT CHARACTERISTICS --- #

# Set for specific patient
has_catheter = 1
has_drain = 1
has_narc_use = 1
duration_IV = 24
post_op_APD = 38
sec_APD = 12
pct_improve_2nd = 0.5

df = tibble(Who_indicated4=0, Surgeon6=0, Return_DietOther=0,
            Catheter1=has_catheter, Drain1=has_drain, 
            Narc_use=has_narc_use, Duration_IV=duration_IV,
            Post_op_APD=post_op_APD, sec_APD=sec_APD,
            percent_improve_2nd=pct_improve_2nd)
stopifnot(all(colnames(df) %in% names(lst_cox$bhat)))

# Normalize
X_row = rownames_to_column(data.frame(x=t(df)),'cn')
dat_norm = tibble(cn=names(lst_cox$mu),bhat=lst_cox$bhat,
                  mu=lst_cox$mu,se=lst_cox$se)
X_row = left_join(X_row,dat_norm,by='cn')
X_row = X_row %>% mutate(x_s = (x-mu)/se)
# eta_patient = summarise(X_row,eta=sum(x_s*bhat))

# --- (3) SURVIVAL DIST --- #

# Align covariance matrix
idx_align = match(X_row$cn,names(lst_cox$bhat))
stopifnot(all(lst_cox$bhat[idx_align] == X_row$bhat))
Sigma = lst_cox$Sigma[idx_align,idx_align]

set.seed(1234)
surv_dist = pm_surv(bhat=X_row$bhat, Sigma = Sigma, Eta=exp(lst_cox$eta),
        Y = lst_cox$y, x=matrix(X_row$x_s,nrow=1),
        nsim = 1000,alpha = 0.05)


# --- (4) PLOT IT --- #
gg_km = ggplot(surv_dist,aes(x=time,y=mu)) +
  theme_bw() + geom_line() +
  labs(y='Survival probability',x='Months',subtitle = 'Shaded area is 95% CI') +
  geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.5) +
  scale_y_continuous(limits=c(0,1))
gg_km
