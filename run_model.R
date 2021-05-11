######################################################
# ---- SCRIPT TO FIT UNIVARIATE SURVIVAL MODELS ---- #

pckgs <- c('magrittr','stringr','dplyr','forcats','tibble',
           'cowplot','ggplot2',
           'survival','flexsurv',
           'glmnet','selectiveInference','mvtnorm')
for (pp in pckgs) { library(pp,character.only = T)}

user <- str_split(getwd(),'\\/')[[1]][3]
stopifnot(user %in% c('erik drysdale','lauren edrman'))
if (user == 'erik drysdale') {
  dir_base <- "C:/Users/erik drysdale/Documents/projects/Pyeloplasty"
} else {
  dir_base <- "C:/Users/lauren erdman/Desktop/pyloplasty"
}
setwd(dir_base)
dir_data <- file.path(dir_base, 'data')

auroc <- function(score,y){
  cls <- y == 1
  n1 <- sum(!cls)
  n2 <- sum(cls)
  U <- sum(rank(score)[!cls])-n1*(n1+1)/2;
  return(1-U/n1/n2);
}

# s <- c(1,3,3)
# l <- c(0,0,1)
# mltools::auc_roc(s,l)
# auroc(s,l)

# 6-12 months different than 1-year

###########################################
# --------- (1) LOAD THE DATA ----------- #

# Find most recent file
udates <- list.files(dir_data) %>% str_subset('\\.rds$') %>% str_split_fixed('\\_',4) %>% extract(,4) %>% as.Date %>% unique
udates <- sort(udates, decreasing = T)
udate <- udates[1]
fn1 <- paste0('pyloplasty_preproc_y_',udate,'.rds')
fn2 <- paste0('pyloplasty_preproc_X_',udate,'.rds')
ydat <- readRDS(file.path(dir_data, fn1)) %>% dplyr::as_tibble() %>% rename(reop=Reoperation,t2e=Time_to_event_allmo)
# Note factors will be expanded for in univariate model
xdat <- readRDS(file.path(dir_data, fn2)) %>% dplyr::as_tibble()
dim(xdat)
# factor lump if less than 2%
cn_fctr <- sapply(xdat,class) %>% data.frame %>% set_colnames('cc') %>% rownames_to_column('cn') %>% 
                filter(cc=='factor') %>% pull(cn)
fprop <- 0.02
xdat <- xdat %>% mutate_at(cn_fctr,function(x) fct_lump(x,prop=fprop))
# Remove NG_tube & IV Vluid
xdat <- xdat %>% dplyr::select(-c(NG_tube,IV_Fluid))
# Impute Blocks
mode_Blocks <- table(xdat$Blocks) %>% sort(T) %>% extract(1) %>% names
xdat$Blocks <- ifelse(is.na(xdat$Blocks),mode_Blocks, as.character(xdat$Blocks))
# Aggregate approach
xdat$Approach <- fct_lump(xdat$Approach,n=2)
# Missing value imputation for APD
cn_apd <- c('Pre_op_APD','Post_op_APD','sec_APD','Last_APD',
            'percent_improve', 'percent_improve_2nd', 'percent_improve_lastffup')
X_apd <- xdat[,cn_apd]
cn_apd <- names(sort(apply(is.na(X_apd),2,sum)))
# Impute lowest missing with median value, then train iterative regression models
for (jj in seq(length(cn_apd))) {
  cn <- cn_apd[jj]
  yy <- pull(X_apd,cn)
  if (jj == 1) {
    print('Median imputation')
    yy <- ifelse(is.na(yy), median(yy,na.rm = T), yy)
    X_apd[,cn] <- yy
  } else {
    print('Parametric imputation')
    ff = formula(str_c(cn_apd[jj],str_c(cn_apd[1:jj-1],collapse='+'),sep='~'))
    mdl <- lm(ff,data=X_apd)
    print(sprintf('Adjusted R-squard: %0.1f%%, DoF: %i',
                  summary(mdl)$adj.r.squared*100,mdl$df.residual))
    yy <- ifelse(is.na(yy), predict(mdl,X_apd), yy)
    X_apd[,cn] <- yy
  }
}
xdat[,cn_apd] <- X_apd[,cn_apd]

################################################
# --------- (2) TRAIN A CURE MODEL ----------- #

# Cured == >30 months
# Not-cured == reop==1
y_cured <- ydat %>% 
  mutate(cured=ifelse(reop==1,'not-cured', ifelse(t2e >= 30, 'cured', 'unknown'))) %>% 
  dplyr::select(c(reop,t2e,cured))
print(table(y_cured$cured))
# one-hot encode x-matrix
Xmat <- model.matrix(~., data=xdat)[,-1]
stopifnot(nrow(Xmat) == nrow(xdat))  # Ensure no missing values
cn_drop <- names(which(apply(Xmat,2,var) < 0.01))
print(sprintf('Removing %i columns for low variance: %s',
              length(cn_drop),str_c(cn_drop,collapse=', ')))
stopifnot(length(cn_drop)==0)
# Subset to y_cured labels
idx_keep <- which(y_cured$cured!='unknown')
# Create glmnet-friendly dataformat
X_cure <- Xmat[idx_keep,]
# Remove low variance
cn_drop_cure = names(which(apply(X_cure,2,var) < 0.01))
sprintf('dropping cure specific columns: %s',cn_drop_cure)
X_cure = X_cure[,setdiff(colnames(X_cure),cn_drop_cure)]
y_cure <- ifelse(y_cured[idx_keep,]$cured=='cured',1,0)
X_cure_s <- scale(X_cure)
mu_X_cure <- attr(X_cure_s,'scaled:center')
se_X_cure <- attr(X_cure_s,'scaled:scale')

# Use CV.glmnet for fit LOO
stime <- Sys.time()
mdl_cv <- cv.glmnet(x=X_cure_s,y=y_cure, family='binomial',
                    nfolds = nrow(X_cure), keep=T,
                    type.measure='deviance', standardize=F, grouped = F)
auroc_cvfold <- apply(mdl_cv$fit.preval, 2, function(eta) auroc(eta,y_cure))
print(Sys.time() - stime)
# Find winning lambda
idx_best <- which.max(auroc_cvfold)
lam_best <- mdl_cv$lambda[idx_best]
# Distribution of AUROC
eta_auroc <- mdl_cv$fit.preval[,which(mdl_cv$lambda==lam_best)]
dist_auroc <- replicate(1000,{
  idx <- sample(length(y_cure),length(y_cure),T)  
  auroc(eta_auroc[idx],y_cure[idx])
})
ggplot(data.frame(x=dist_auroc),aes(x=x)) + geom_histogram(fill='grey',color='red') + 
  geom_vline(xintercept = quantile(dist_auroc,0.025)) + 
  geom_vline(xintercept = quantile(dist_auroc,0.975)) + 
  theme_bw() + ggtitle('Bootstrap distribution of LOO-AUROC for Cure Model')
print(max(auroc_cvfold))
quantile(dist_auroc,c(0.025,0.975))

# Refit lasso with full data
mdl_cure <- glmnet(x=X_cure_s,y=y_cure, family='binomial',
                   lambda = lam_best, standardize = F)

# Run selective inference
SI_cure <- fixedLassoInf(x=X_cure_s,y=y_cure,
                         family='binomial',alpha=0.05,
                         beta = as.vector(coef(mdl_cure)),
                         lambda = lam_best*nrow(X_cure))
bhat_cure <- tibble(cn=names(SI_cure$vars),coef=SI_cure$coef0,
       pval=SI_cure$pv,lb=SI_cure$ci[,1],ub=SI_cure$ci[,2])
bhat_cure %>% filter(pval < 0.05)
# Scale Xmat using cure params
Xmat_scure <- sweep(sweep(Xmat[,colnames(X_cure_s)],2,mu_X_cure,'-'),2,se_X_cure,'/')

# Get the predicted weights
cure_weights <- 1/(1+exp(-as.vector(predict(mdl_cure, newx=Xmat_scure))))

y_cured <- y_cured %>% mutate(cweights=cure_weights) %>% 
  mutate(cweights2=ifelse(cured=='not-cured',0,cweights))

##############################################
# --------- (2) FIT HIGH-DIM COX ----------- #

p <- seq(0.5,1,0.01)
dat_p <- tibble(p=p,m=sapply(p, function(x) mean(y_cured$cweights<x)))
dat_p %>% mutate(dd=m-lag(m,1)) %>% tail(10)
plot(dat_p$p, dat_p$m)
abline(v=0.955)

idx_surv <- which((y_cured$cweights2 < 0.95) & (y_cured$cured!='cured'))
print(sprintf('Using %i of %i non-cured rows', length(idx_surv),nrow(y_cured)))

y_cured %>% mutate(is_w = cweights2 < 0.95) %>% 
  filter(cured != 'cured') %>% group_by(cured,is_w) %>% count

# Remove patients who we know to be cured
X_surv <- Xmat[idx_surv,]
# Remove low variance
cn_drop_surv = names(which(apply(X_surv,2,var) < 0.01))
sprintf('dropping surv specific columns: %s',cn_drop_surv)
X_surv = X_surv[,setdiff(colnames(X_surv),cn_drop_surv)]
X_surv_s <- scale(X_surv)
mu_X_surv <- attr(X_surv_s,'scaled:center')
se_X_surv <- attr(X_surv_s,'scaled:scale')
y_surv <- with(y_cured[idx_surv,],Surv(t2e, reop))

# Use CV.glmnet for fit LOO
stime <- Sys.time()
cv_surv <- cv.glmnet(x=X_surv_s,y=y_surv, family='cox',
                    nfolds = nrow(X_surv), keep=T,
                    type.measure='deviance', standardize=F)
print(Sys.time() - stime)
res_conc <- apply(cv_surv$fit.preval, 2, 
                  function(eta) survConcordance.fit(y=y_surv, x=eta)) %>% 
  t %>% as_tibble %>% mutate(num=concordant+0.5*tied.risk) %>% 
  mutate(den=num+discordant) %>% mutate(conc=num/den,lam=cv_surv$lambda) %>% 
  dplyr::select(c(lam,conc))
lam_surv_star <- res_conc %>% arrange(-conc) %>% pull(lam) %>% head(1)
eta_star <- cv_surv$fit.preval[,which(cv_surv$lambda==lam_surv_star)]

dist_conc <- replicate(1000,{
  idx <- sample(nrow(y_surv),nrow(y_surv),T)  
  survConcordance(y_surv[idx]~ eta_star[idx])$concordance
})
ggplot(data.frame(x=dist_conc),aes(x=x)) + geom_histogram(fill='grey',color='blue') + 
  geom_vline(xintercept = quantile(dist_conc,0.025)) + 
  geom_vline(xintercept = quantile(dist_conc,0.975)) + 
  theme_bw() + ggtitle('Bootstrap distribution of LOO-Concordance for Survival Model')

res_conc %>% arrange(-conc) %>% head(1) %>% print
quantile(dist_conc,c(0.025,0.975))


# Refit
mdl_surv <- glmnet(x=X_surv_s,y=y_surv, family='cox', standardize=F,
                   lambda = lam_surv_star)
# Get SI
SI_surv <- fixedLassoInf(x=X_surv_s,y=y_surv[,1],status = y_surv[,2],
                         family='cox',alpha=0.05,
                         beta = as.vector(coef(mdl_surv)),
                         lambda = lam_surv_star*nrow(X_surv))
bhat_surv <- tibble(cn=colnames(X_surv_s)[SI_surv$vars] ,coef=SI_surv$coef0,
                    pval=SI_surv$pv,lb=SI_surv$ci[,1],ub=SI_surv$ci[,2])
bhat_surv %>% filter(pval < 0.1)

############################################
# --------- (3) FIT GLM MODELS ----------- #

df_bhat <- rbind(mutate(bhat_cure,tt='cure'),mutate(bhat_surv,tt='cox'))
df_bhat <- df_bhat %>% 
  mutate(is_sig=ifelse(pval<0.05,T,F)) %>% 
  mutate_at(vars(c('lb','ub')),list(~ifelse(abs(.)==Inf,NA, .))) %>% 
  mutate(bound=ifelse(sign(coef)==1, lb, ub)) %>%
  mutate(bound=ifelse(is_sig, bound, NA))

X_surv_s_sub = as.data.frame(X_surv_s[,bhat_surv$cn])
glm_surv = coxph(y_surv~.,data=X_surv_s_sub)
glm_cure = glm(y_cure~X_cure_s[,bhat_cure$cn],family='binomial')
glm_bhat = rbind(tibble(tt='cox',mdl='glm',cn=names(coef(glm_surv)),coef=coef(glm_surv)),
                 tibble(tt='cure',mdl='glm',cn=bhat_cure$cn,coef=coef(glm_cure)[2:length(coef(glm_cure))]))
tmp1 = mdl_surv %>%
  coef %>% as.matrix %>% extract(,1) %>% data.frame(coef=.) %>% 
  rownames_to_column('cn') %>% as_tibble() %>% mutate(tt='cox')
tmp2 = mdl_cure %>% coef %>% as.matrix %>% extract(,1) %>% 
  data.frame(coef=.) %>% 
  rownames_to_column('cn') %>% as_tibble() %>% 
  mutate(tt='cure') %>% filter(cn != '(Intercept)')
lasso_bhat = mutate(rbind(tmp1, tmp2),mdl='lasso')

all_bhat = rbind(glm_bhat, lasso_bhat,mutate(df_bhat[,c('tt','cn','coef')],mdl='SI'))
all_bhat = filter(all_bhat, coef != 0)
all_bhat %>% filter(tt=='cox') %>% tidyr::pivot_wider(cn,names_from='mdl',values_from = 'coef')

#########################################
# --------- (4) PM EXAMPLES ----------- #

coef_cox = as.vector(glm_surv$coefficients)
risk_cox = exp(as.matrix(X_surv_s_sub) %*% coef_cox)
id_mi = which.min((X_surv_s_sub$Post_op_APD-1)**2)
id_mx = which.min((X_surv_s_sub$Post_op_APD+1)**2)

X_mi = X_surv_s_sub[id_mi,]
X_mx = X_surv_s_sub[id_mx,]
risk_mi = risk_cox[id_mi,]
risk_mx = risk_cox[id_mx,]
y_mi = y_surv[id_mi,,drop=F]
y_mx = y_surv[id_mx,,drop=F]

# FUNCTION WRAPPER TO GET SURV CIs from SIMULATION
pm_surv = function(mu, Sigma, X, Y, x, y, nsim=100, alpha=0.05) {
  set.seed(nsim)
  if (!is.matrix(X)) { X = as.matrix(X) }
  if (!is.matrix(x)) { x = as.matrix(x) }
  mu_sim = rmvnorm(n=nsim, mean=mu, sigma=Sigma)
  holder = list()
  Eta = exp(X %*% mu)
  for (i in seq(nsim)) {
    mu_i = mu_sim[i,]
    eta_i = exp(x %*% mu_i)[1,1]
    res_i = coxsurv.fit(ctype=1, stype = 1, se.fit = FALSE, cluster = NULL,
                varmat = Sigma, y = Y, x = X, wt = rep(1,nrow(X)),
                risk = Eta, y2 = y, x2 = x, risk2 = eta_i,
                position = NULL, strata = NULL, oldid = NULL,
                strata2 = NULL, id2 = NULL,unlist = TRUE)
    res_i = mutate(as.data.frame(do.call('cbind',res_i[c('time','surv')])),idx=i)
    holder[[i]] = res_i
  }
  res_sim = as_tibble(do.call('rbind',holder))
  res_sim = res_sim %>% group_by(time) %>%  
    summarise(mu=mean(surv),lb=quantile(surv,alpha/2),ub=quantile(surv,1-alpha/2)) %>%
    arrange(time)
  return(res_sim)
}
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Individualized curves
alpha = 0.05
tmp1 = pm_surv(mu=coef_cox, Sigma=glm_surv$var, X=X_surv_s_sub, Y=y_surv,
        x=X_mi,y=y_mi,nsim=1000,alpha=alpha)
tmp2 = pm_surv(mu=coef_cox, Sigma=glm_surv$var, X=X_surv_s_sub, Y=y_surv,
               x=X_mx,y=y_mx,nsim=1000,alpha=alpha)
zscore = qnorm(1-alpha/2)
tmp3 = survfit(y_surv~1,conf.int=1-alpha, se.fit=T)
tmp3 = tibble(time=tmp3$time,mu=tmp3$surv,se=tmp3$std.err,tt='KM')
tmp3 = tmp3 %>% mutate(lb=mu-zscore*se, ub=mu-zscore*se) %>% dplyr::select(-se)
df_km = rbind(mutate(tmp1,tt='Low'),mutate(tmp2,tt='High'),tmp3)

colz = c(gg_color_hue(2)[1],'black',gg_color_hue(2)[2])
gg_km = ggplot(filter(df_km,tt!='KM'),aes(x=time,y=mu,color=tt,fill=tt)) + 
  theme_bw() + geom_line() + 
  labs(y='Survival probability',x='Months',subtitle = 'Shaded area is 95% CI') + 
  scale_color_discrete(name='APD') +
  scale_fill_discrete(name='APD') +
  geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.5) + 
  scale_y_continuous(limits=c(0,1))
gg_km

########################################
# --------- (5) MAKE PLOTS ----------- #

posd = position_dodge(0.5)
# SI==GLM for Cox, different for cure
gg_comp = ggplot(all_bhat,aes(x=cn,y=coef,color=mdl)) + theme_bw() + 
  geom_point(position = posd) + 
  labs(y='Coefficient') + geom_hline(yintercept = 0) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=90)) + 
  facet_wrap(~tt)
gg_comp

gg_bhat <- ggplot(df_bhat, aes(x=fct_reorder2(cn,coef,is_sig), y=coef,color=is_sig)) + 
  geom_point(size=3) + facet_wrap(~tt,scales='free_x') + 
  theme_bw() + ggtitle('Significance for SI coefficients') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=90)) + 
  geom_hline(yintercept = 0, linetype='dashed')
  # geom_linerange(aes(ymin=bound,ymax=coef)) + 
  # geom_linerange(aes(ymin=coef,ymax=bound))
gg_bhat

#######################################################
# --------- (5) SAVE FIGURES FOR NOTEBOOK ----------- #

df_auroc <- data.frame(auroc=dist_auroc)
df_conc <- data.frame(conc=as.vector(dist_conc))
df_bhat <- df_bhat %>% dplyr::select(-c(lb,ub,bound))
save(dat_p,df_auroc, df_conc, df_bhat, df_km, file = file.path(dir_base,'fig_data.RData'))
     

