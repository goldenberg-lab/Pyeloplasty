######################################################
# ---- SCRIPT TO FIT UNIVARIATE SURVIVAL MODELS ---- #

pckgs <- c('magrittr','stringr','dplyr','forcats','tibble',
           'cowplot','ggplot2',
           'survival','flexsurv',
           'glmnet')
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
  n1<-sum(!cls); sum(cls)->n2;
  U<-sum(rank(score)[!cls])-n1*(n1+1)/2;
  return(1-U/n1/n2);
}

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
# factor lump if less than 2%
cn_fctr <- sapply(xdat,class) %>% data.frame %>% set_colnames('cc') %>% rownames_to_column('cn') %>% 
                filter(cc=='factor') %>% pull(cn)
fprop <- 0.02
xdat <- xdat %>% mutate_at(cn_fctr,function(x) fct_lump(x,prop=fprop))
# Remove NG_tube & IV Vluid
xdat <- xdat %>% select(-c(NG_tube,IV_Fluid))
# Impute Blocks
mode_Blocks <- table(xdat$Blocks) %>% sort(T) %>% extract(1) %>% names
xdat$Blocks <- ifelse(is.na(xdat$Blocks),mode_Blocks, as.character(xdat$Blocks))
# Aggregate approach
xdat$Approach <- fct_lump(xdat$Approach,n=2)

# --------- (2) TRAIN A CURE MODEL ----------- #

# Cured == >30 months
# Not-cured == reop==1
y_cured <- ydat %>% mutate(cured=ifelse(reop==1,'not-cured', ifelse(t2e >= 30, 'cured', 'unknown'))) %>% 
                select(c(reop,t2e,cured))
# one-hot encode x-matrix
Xmat <- model.matrix(~., data=xdat)[,-1]
# Xmat <- cbind(t2e=y_cured$t2e,Xmat)
cn_drop <- names(which(apply(Xmat,2,var) < 0.01))
stopifnot(length(cn_drop)==0)
# Subset to y_cured labels
idx_keep <- which(y_cured$cured!='unknown')

y_cured[idx_keep,] %>% head
Xmat[idx_keep,] %>% head(1) %>% t

jj <- 0
phat <- rep(NA, length(idx_keep))
for (ii in idx_keep) {
  jj <- jj + 1
  if (jj %% 50 == 0) {
    print(jj)
  }
  itrain <- idx_keep[-jj]
  iX <- Xmat[itrain,]
  mu_X <- apply(iX,2,mean)
  se_X <- apply(iX,2,sd)
  iX <- sweep(sweep(iX,2,mu_X,'-'),2,se_X,'/')
  iy <- ifelse(y_cured[itrain,]$cured=='cured',1,0)
  # Fit logistic lasso
  ilasso <- glmnet(x=iX,y=iy,family='binomial',nlambda=100,standardize=F)
  eta_mat <- as.matrix(cbind(1,iX) %*% coef(ilasso))
  phat_mat <- 1/(1+exp(-eta_mat))
  res_mat <- iy - phat_mat
  e2 <- apply(res_mat, 2, function(x) sqrt(sum(x**2)))
  z2 <- nrow(iX)*ilasso$lambda / e2
  thresh <- qnorm(1-0.05/ncol(iX))
  istar <- which.min((z2 - thresh)^2)
  ibhat <- coef(ilasso)[,istar]
  phat[jj] <- as.numeric(c(Intercept=1,(Xmat[ii,]-mu_X)/se_X) %*% ibhat)
}

auc_roc(phat,ifelse(y_cured[idx_keep,]$cured=='cured',1,0))

X_cure <- Xmat[idx_keep,]
y_cure <- ifelse(y_cured[idx_keep,]$cured=='cured',1,0)

set.seed(1234)
mdl_cv <- cv.glmnet(x=X_cure,y=y_cure, family='binomial',
                    nfolds = 5, keep=T,type.measure='auc')
auroc_cvfold <- apply(mdl_cv$fit.preval, 1, function(eta) auc_roc(eta,y_cure))
max(auroc_cvfold)

# Get the predicted weights
cure_weights <- cbind(1,Xmat) %*% coef(mdl_cv,s='lambda.min')
cure_weights <- 1/(1+exp(-cure_weights[,1]))
rm(list=c('X_cure','y_cure'))

y_cured <- y_cured %>% mutate(cweights=cure_weights)
y_cured <- y_cured %>% mutate(weights=ifelse(cured=='cured',0,ifelse(cured=='not-cured',1,1-cweights)))

# --------- (2) FIT UNIVARIATE EXPONENTIAL MODELS ----------- #

# Remove patients who we know to be cured
X_sub <- xdat[y_cured$weights>0,]
cn_numeric <- names(which(sapply(X_sub,class)=='numeric'))
# Scale the continuous variables to help with likelihood
X_sub[,cn_numeric] <- scale(X_sub[,cn_numeric])

y_sub <- y_cured %>% filter(weights >0) %>% select(-cured)
So_sub <- with(y_sub, Surv(t2e, reop))

cn_X <- colnames(X_sub)

# Gets SEs close to zero
# ff <- as.formula(str_c('So_sub~',str_c(colnames(X_sub),collapse='+')))
# flexsurvreg(ff, data=X_sub, dist='exponential',weights=y_sub$weights)

holder <- vector('list',length(cn_X))
jj <- 0
for (cn in cn_X) {
  print(cn)
  jj <- jj + 1
  # Exponential model
  fcn <- as.formula(paste0('So_sub~',cn))
  mdl_exp <- flexsurvreg(fcn, data=X_sub, dist='exponential',weights=y_sub$weights)
  # Get the p-values
  vv <- names(mdl_exp$coefficients)
  tmp <- as_tibble(mdl_exp$res) %>% mutate(cn=cn,vv=vv,pval=2*pnorm(abs(est/se),lower.tail=F)) %>% 
    select(c(cn,vv,est,pval))
  holder[[jj]] <- tmp
}
dat_exp <- do.call('rbind',holder)

# Plot the distribution of 
dat_exp %>% filter(vv!='rate') %>% arrange(pval) %>% mutate(pfdr=p.adjust(pval,'bonferroni'))

# bhat_exp <- mdl_exp$coefficients
# lam_exp <- exp(X[ii,,drop=F] %*% bhat_exp)
# mat[ii,'med_exp']  <- qexp(p=0.5,rate=lam_exp)

# --------- (3) COMPARE TO SURVIVAL MODELS ----------- #


So_all <- with(ydat,Surv(t2e,reop))
X_all <- model.matrix(~.,data=xdat)[,-1]
set.seed(1234)
cox_all <- cv.glmnet(x=X_all, y=So_all, type.measure='C', nfolds=5, family='cox', keep=T)
conc_all <- apply(cox_all$fit.preval, 2, function(x) concordancefit(y=So_all,x)$concordance)
max(conc_all)


X_cure <- model.matrix(~.,data=X_sub)[,-1]
y_cure <- So_sub
set.seed(1234)
cox_cure <- cv.glmnet(x=X_cure, y=y_cure, type.measure='C', nfolds=5, family='cox', keep=T, weights=y_sub$weights)
conc_cure <- apply(cox_cure$fit.preval, 2, function(x) concordancefit(y=y_cure,x)$concordance)
max(conc_cure)


# --------- (4) PLOTS ----------- #

# Notice the survival flatline
plot(survfit(Surv(t2e, reop) ~ 1, data = ydat),
     xlab = "Months", ylab = "Overall survival probability",main='KM curve with no weights')

plot(survfit(Surv(t2e, reop) ~ 1, data = y_cured, weights = weights, subset = weights>0),
     xlab = "Months", ylab = "Overall survival probability",main='KM Curve with Weights')

plot(survfit(Surv(t2e, reop) ~ 1, data = y_cured, subset = reop==1),
     xlab = "Months", ylab = "Overall survival probability",main='KM Curve for event-only')

gg_weights <- ggplot(y_cured,aes(x=cweights)) + theme_bw() +
  geom_histogram(bins=10,color='red',fill='grey') + facet_wrap(~cured) + 
  ggtitle('Distribution of cure probabilities')
gg_weights

gg_rate <- ggplot(filter(dat_exp,vv=='rate' & pval<0.05),aes(x=est)) + theme_bw() + 
            geom_histogram(bins=10,color='blue',fill='grey') + 
            ggtitle('Distribution of univariate rate parameters')
gg_rate

