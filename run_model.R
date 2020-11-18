#############################################
# ---- SCRIPT TO LOAD DATA AND ANALYZE ---- #

library(survival)
library(flexsurv)
library(ggplot2)
library(cowplot)
library(dplyr)
library(magrittr)
library(readxl)

dir_base <- "C:/Users/erik drysdale/Documents/projects/Pyeloplasty"
setwd(dir_base)

# "Blocks"
# "Surgeon"
pred_cols = c('Reoperation','Time_to_event_allmo','Definite_success',
  "Sex","Age_sx_Mos","Who_indicated","Sex_Provider",  
              "Side","Approach","Intraop_finding",
              "Angiopexy","Salle_Stent", "JJ_Stent",
              "OR_Time","Nx_Vx_Prophy", "Dex" ,
              "NG_tube","Catheter" , "Drain"   ,
              "Narc_Rx","Narc_use",
              "PCA_use","Epidural",
              "ketoralac", "Duration_IV","Constipation_prophy","Mobilization") 

data <- readRDS('pyloplasty_preproc.rds') %>% dplyr::as_tibble() %>% select(pred_cols)
df_sub <- data %>% rename(t2e='Time_to_event_allmo',reop='Reoperation',
                age='Age_sx_Mos',cured='Definite_success') %>% 
  mutate(cured=ifelse(is.na(cured),0,1))

# Cure is 30 months! Drop reop if >30 months
# Exclude Angiopexy



# Columns to use
# cn_use <- c('Reoperation','Time_to_event_allmo','Age_sx_Mos','Definite_success')

# df <- read_xlsx(path='Pyeloplasty.xlsx')
# df_sub <- df %>% select(cn_use) %>% rename(t2e='Time_to_event_allmo',reop='Reoperation',
#                                            age='Age_sx_Mos',cured='Definite_success') %>% 
#   mutate(cured=ifelse(is.na(cured),0,1))
# df_sub <- df_sub %>% filter(cured==0)

# mdl_exp <- flexsurvreg(Surv(t2e,reop)~1,data=df_sub,dist='exponential')
# qexp(p=0.5,rate=exp(mdl_exp$coefficients)) / 12
# 
# mdl_wei <- flexsurvreg(Surv(t2e,reop)~1,data=df_sub,dist='weibull')
# bhat_wei <- mdl_wei$coefficients
# shape_trans <- exp(bhat_wei['shape'])
# scale_trans <- (1/exp(bhat_wei['scale']))^shape_trans
# qweibullPH(p=0.5, shape=shape_trans, scale=scale_trans) / 12
# 
# plot(survfit(Surv(t2e, reop) ~ 1, data = df_sub), 
#      xlab = "Months", 
#      ylab = "Overall survival probability")


###########################################################
# ---- Step 1: Do leave one-out for cross validation ---- #

X <- model.matrix(~.,data=select(df_sub,-c(t2e,reop,cured)))

So <- with(df_sub,Surv(t2e,reop))
n <- nrow(df_sub)

cn_mat <- c('med_exp', 'med_wei')
mat <- matrix(nrow=n, ncol=length(cn_mat), dimnames = list(NULL, cn_mat))

for (ii in seq(n)) {
  # Exponential model
  mdl_exp <- flexsurvreg(So[-ii]~X[-ii,-1],dist='exponential')
  bhat_exp <- mdl_exp$coefficients
  lam_exp <- exp(X[ii,,drop=F] %*% bhat_exp)
  mat[ii,'med_exp']  <- qexp(p=0.5,rate=lam_exp)
  
  # # Weibull model
  # mdl_wei <- flexsurvreg(So[-ii]~X[-ii,-1],dist='weibull')
  # bhat_wei <- mdl_wei$coefficients
  # shape_trans <- exp(bhat_wei['shape'])
  # bvec_trans <- log((1/exp(bhat_wei[2:length(bhat_wei)]))^shape_trans)
  # scale_trans <- exp(as.vector(X[ii,,drop=F] %*% bvec_trans))
  # mat[ii,'med_wei'] <- qweibullPH(p=0.5, shape=shape_trans, scale=scale_trans)
    
}

nsim <- 250 
bs_hold <- rep(NA, nsim)
for (jj in seq(nsim)) {
  idx <- sample(seq(n),replace = T)
  score_jj <- mat[idx,]
  So_jj <- So[idx]
  res <- apply(score_jj, 2, function(z) survival:::concordancefit(y=So_jj, x=-z)$concordance)
  bs_hold[jj] <- res[1] - res[2]
}
hist(bs_hold)



