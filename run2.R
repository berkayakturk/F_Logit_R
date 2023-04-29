# load the auxilary functions

source("functions.R")
source("dgp3.R")
source("dgp4.R")

ccr_pls = numeric() 
ccr_pca = numeric()
ccr_rpls = numeric()
ccr_rpca = numeric() 

auc_pls = numeric() 
auc_pca = numeric()
auc_rpls = numeric() 
auc_rpca = numeric()

bhat_pls = list()
bhat_pca = list()
bhat_rpls = list()
bhat_rpca = list() 

rangeval = c(1,21) # # value range
gpt = seq(1, 21, length.out = 101) 
out_p = 0.1 # outlier proportion 10%
nsim = 10

# data generation

set.seed(123)
dat = dgp3(350, 350)
X = dat$X
Y = c(rep(0,350), rep(1,350))
dat_test = dgp3(150, 150)
X_test=dat_test$X
Y_test=c(rep(0,150), rep(1,150))

# outlier generation

dat2 = dgp4(350, 350)
X_out = dat2$X
    
index = which(Y == 1)
out_index=sample(index, out_p * 700, replace = F)
    
X[out_index,] = X_out[out_index,]
Y[out_index] = 0

# get the desing matrix for fpls (H)

H = getHmat(data = X, nbasis = 15, rangeval = rangeval)
Hmat=H$Hmat
Bsf=H$Bsf
    
H2 = getHmat(data = X_test, nbasis = 15, rangeval = rangeval)
Hmat_test = H2$Hmat

# run the pls, pca, and rpls and rpca methods
    
pls_model = log_fpls2(Y = Y, Hmat = Hmat, Hmat_test = Hmat_test, alpha = 0.05, hmax = 10, Bsf=Bsf,  model = "ls")
pca_model = log_rpca(Y = Y, X = X, X_test = X_test, nbasis = 15, ncomp = 5, gp = gpt, emodel = "classical")
rpls_model = log_fpls2(Y = Y, Hmat = Hmat, Hmat_test = Hmat_test, alpha = 0.05, hmax = 10, Bsf=Bsf,  model = "wle")
rpca_model = log_rpca(Y = Y, X = X, X_test = X_test, nbasis = 15, ncomp = 5, gp = gpt, emodel = "robust")

# nbasis = The number of b spline basis functions
# ncomp = Number of principal components
# gp = A Gaussian process range
# emodel or model = A character string specifying the estimation model to be used.
# alpha = The significance level for variable selection, e.g., 0.05
# hmax = The maximum number of steps for the FPLS algorithm.
# Bsf = The matrix of basis functions used for expanding the variables.
    
# ccr values

pls_ccr = ccr(Y_test, pls_model$preds)
pca_ccr = ccr(Y_test, pca_model$preds)
rpls_ccr = ccr(Y_test, rpls_model$preds)
rpca_ccr = ccr(Y_test, rpca_model$preds)

ccr_pls = pls_ccr$ccr # 0.68
ccr_pca = pca_ccr$ccr # 0.9
ccr_rpls = rpls_ccr$ccr # 0.93
ccr_rpca = rpca_ccr$ccr # 0.9266667

# auc values

auc_pls = auc(Y_test ~ pls_model$preds)[1] # 0.68
auc_pca = auc(Y_test ~ c(pca_model$preds))[1] # 0.9
auc_rpls = auc(Y_test ~ rpls_model$preds)[1] # 0.93
auc_rpca = auc(Y_test ~ c(rpca_model$preds))[1] # 0.9266667

# predicted parameter functions
    
bhat_pls = c(pls_model$beta_hat_s)
bhat_pca = c(pca_model$beta_hat_s)
bhat_rpls = c(rpls_model$beta_hat_s)
bhat_rpca = c(rpca_model$beta_hat_s)    
 
    