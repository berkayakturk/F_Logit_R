#___________Packages__________

library(fda) # functional data analysis
library(robustX) # Robust estimation and inference for regression and scale  
library(matrixStats) # Functions that Apply to Rows and Columns of Matrices
library(pROC) # Display and analyze ROC curves
library(pcaPP) # Robust PCA methods 
library(robustbase) # Basic Robust Statistics 
library(RobStatTM) # Robust Statistics: Bianco and Yohai estimator for logistic regression
library(expm) # Matrix exponential, logarithm, square root, and related functions
library(wle) # Weighted Likelihood Estimation

#____________________________

options(warn = -1)


#__________________________________Logistics functional principal component analysis__________________________________

getPCA <-
  function(data, nbasis, ncomp, gp, emodel = c("classical", "robust")){
    #pcaPP paketini kurman gerek
    emodel <- match.arg(emodel)
    n <- dim(data)[1]
    p <- dim(data)[2]
    dimnames(data) = list(as.character(1:n), as.character(1:p))
    bs_basis <- create.bspline.basis(rangeval = c(gp[1], gp[p]), nbasis = nbasis)
    inp_mat <- inprod(bs_basis, bs_basis)
    sinp_mat <- sqrtm(inp_mat)
    evalbase = eval.basis(gp, bs_basis)
    fdobj <- fdPar(bs_basis, int2Lfd(2), lambda=0)
    pcaobj <- smooth.basisPar(gp, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
    
    if(emodel == "classical"){
      
      mean_coef <- apply(t(pcaobj$coefs), 2, mean)
      sdata <- scale(t(pcaobj$coefs), scale = FALSE)
      new.data <- sdata %*% sinp_mat
      dcov <- cov(new.data)
      d.eigen <- eigen(dcov)
      loads <- d.eigen$vectors[,1:ncomp]
      PCs <- solve(sinp_mat) %*% loads
      colnames(PCs) = 1:ncomp
      for(i in 1:ncomp)
        colnames(PCs)[i] = paste("PC", i, sep = "")
      PCAcoef <- fd(PCs, bs_basis)
      mean_coef <- fd(as.vector(mean_coef), bs_basis)
      pcaobj2 <- pcaobj
      pcaobj2$coefs <- t(sdata)
      PCAscore <- inprod(pcaobj2, PCAcoef)
      colnames(PCAscore) = 1:ncomp
      for(i in 1:ncomp)
        colnames(PCAscore)[i] = paste("Score", i, sep = "")
    }else if(emodel == "robust"){
      mean_coef <- pcaPP::l1median(t(pcaobj$coefs), trace = -1)
      sdata <- scale(t(pcaobj$coefs), center = mean_coef, scale = FALSE)
      new.data <- sdata %*% sinp_mat
      ppur <- PCAproj(new.data, ncomp)
      loads <- ppur$loadings
      PCs <- solve(sinp_mat) %*% loads
      colnames(PCs) = 1:ncomp
      for(i in 1:ncomp)
        colnames(PCs)[i] = paste("PC", i, sep = "")
      PCAcoef <- fd(PCs, bs_basis)
      mean_coef <- fd(as.vector(mean_coef), bs_basis)
      pcaobj2 <- pcaobj
      pcaobj2$coefs <- t(sdata)
      PCAscore <- inprod(pcaobj2, PCAcoef)
      colnames(PCAscore) = 1:ncomp
      for(i in 1:ncomp)
        colnames(PCAscore)[i] = paste("Score", i, sep = "")
    }
    return(list(PCAcoef = PCAcoef, PCAscore = PCAscore, meanScore = mean_coef,
                bs_basis = bs_basis, evalbase = evalbase))
  }



getPCA.test <-
  function(data, bs_basis, PCAcoef, gp, emodel = c("classical", "robust")){
    
    emodel <- match.arg(emodel)
    n <- dim(data)[1]
    p <- dim(data)[2]
    dimnames(data) = list(as.character(1:n), as.character(1:p))
    pcaobj <- smooth.basisPar(gp, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
    
    if(emodel == "classical"){
      sdata <- scale(t(pcaobj$coefs), scale = FALSE)
      pcaobj2 <- pcaobj
      pcaobj2$coefs <- t(sdata)
      PCAscore_test = inprod(pcaobj2, PCAcoef)
      colnames(PCAscore_test) = 1:dim(PCAcoef$coefs)[2]
      for(i in 1:dim(PCAcoef$coefs)[2])
        colnames(PCAscore_test)[i] = paste("Score", i, sep = "")
    }else if(emodel == "robust"){
      mean_coef <- pcaPP::l1median(t(pcaobj$coefs), trace = -1)
      sdata <- scale(t(pcaobj$coefs), center = mean_coef, scale = FALSE)
      pcaobj2 <- pcaobj
      pcaobj2$coefs <- t(sdata)
      PCAscore_test = inprod(pcaobj2, PCAcoef)
      colnames(PCAscore_test) = 1:dim(PCAcoef$coefs)[2]
      for(i in 1:dim(PCAcoef$coefs)[2])
        colnames(PCAscore_test)[i] = paste("Score", i, sep = "")
    }
    
    return(PCAscore.test = PCAscore_test)
  }
#________________________________________________________

#__________________________________Robust FPCA and classical FPCA based on logregBY estimator (2002)__________________________________

log_rpca = function(Y, X, X_test, nbasis, ncomp, gp, emodel){
  
  if(emodel == "classical"){
    
    cpca <- getPCA(data=X, nbasis=nbasis, ncomp=ncomp, gp=gp, emodel = "classical") 
    
    cpca_s <- cpca$PCAscore
    c_bs_basis <- cpca$bs_basis
    c_PCAcoef <- cpca$PCAcoef
    
    cpca_test <- getPCA.test(X_test, c_bs_basis, c_PCAcoef, gpt, emodel = "classical")
    
    c_log_model = glm(Y~cpca_s, family = binomial)
    c_Gam_par = c_log_model$coefficients[-1]
    c_fits = round(c_log_model$fitted.values)
    c_Gam_constant = c_log_model$coefficients[1]
    c_beta_hat = c_PCAcoef$coefs%*%c_Gam_par
    c_beta_hat_s = cpca$evalbase%*%c_beta_hat
    
    c_preds = cbind(1, cpca_test) %*% c(c_Gam_constant, c_Gam_par) 
    c_preds_y = round(exp(c_preds) / (1 + exp(c_preds)))
    
    return(list(fits = c_fits, beta_hat_s = c_beta_hat_s, beta_hat=c_beta_hat,
                Gam_constant=c_Gam_constant, preds=c_preds_y))
    
  }else if(emodel == "robust"){
    
    rpca <- getPCA(data=X, nbasis=nbasis, ncomp=ncomp, gp=gp, emodel = "robust")
    
    rpca_s <- rpca$PCAscore
    r_bs_basis <- rpca$bs_basis
    r_pcacoef <- rpca$PCAcoef
    
    rpca_test <- getPCA.test(X_test, r_bs_basis, r_pcacoef, gpt, emodel = "robust")
    
    r_log_model = logregBY(rpca_s, Y, intercept=1) #1 veya 0, bir engelin dahil edilip edilemeyeceğini gösterir
    r_Gam_par = r_log_model$coefficients[-1]
    r_fits = round(r_log_model$fitted.values)
    r_Gam_constant = r_log_model$coefficients[1]
    
    r_beta_hat = r_pcacoef$coefs%*%r_Gam_par
    r_beta_hat_s = rpca$evalbase%*%r_beta_hat
    
    r_preds = cbind(1, rpca_test) %*% c(r_Gam_constant, r_Gam_par) 
    r_preds_y = round(exp(r_preds) / (1 + exp(r_preds)))
    
    return(list(fits = r_fits, beta_hat_s = r_beta_hat_s, beta_hat=r_beta_hat,
                Gam_constant=r_Gam_constant, preds=r_preds_y))
  } 
}

#________________________________________________________

#____________________L1 median function_________________

med_L1 = function(X){
  nc = dim(X)[2]
  nr = dim(X)[1]
  
  medin = apply(X, 2, median, na.rm = TRUE)
  loc = function(X, nc) X - rep(nc, each = nr)
  scl = apply(abs(loc(X, medin)), 2, mean, trim=0.2, na.rm = TRUE)
  ps = scl>0
  
  if(sum(ps)){
    medin[ps] = L1median(X[, ps, drop=FALSE], m.init = medin[ps],
                         pscale=scl[ps])$estimate
  }
  return(medin)
}
#________________________________________________________


#______________________Obtain H matrix___________________
getHmat = function(data, nbasis, rangeval){
  n = dim(data)[1]
  p = dim(data)[2]
  dimnames(data)=list(as.character(1:n), as.character(1:p))
  grid_points = seq(rangeval[1], rangeval[2], length.out = p)
  bs_basis = create.bspline.basis(rangeval, nbasis = nbasis)
  evalbase = eval.basis(grid_points, bs_basis)
  innp_mat = inprod(bs_basis, bs_basis)
  fdobj = fdPar(bs_basis, int2Lfd(2), lambda=0)
  pcaobj = smooth.basisPar(grid_points, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
  H = t(pcaobj$coefs) %*% innp_mat
  return(list(Hmat = H, Bsf = evalbase))
}
#________________________________________________________


#____________________________Logistic Functional PLS proposed based on Dodge and Whittaker (2009)_____________________________
log_fpls2 = function(Y, Hmat, Hmat_test, alpha, hmax, Bsf, model){
  
  Htest = Hmat_test
  
  if(model == "ls"){
    Hmat = scale(Hmat)
    Hmat_test = scale(Hmat_test)
  }else  if(model == "wle"){
    minw = 0.1
    medH = med_L1(Hmat)
    medH_test = med_L1(Hmat_test)
    cH = sapply(1:ncol(Hmat), function(i) (Hmat[,i]-medH[i]))
    cH_test = sapply(1:ncol(Hmat_test), function(i) (Hmat_test[,i]-medH_test[i]))
    Hmat = sapply(1:ncol(cH), function(i) cH[,i]/median(abs(cH[,i]))/0.6745)
    Hmat_test = sapply(1:ncol(cH_test), function(i) cH_test[,i]/median(abs(cH_test[,i]))/0.6745)
    wH = sapply(1:ncol(Hmat), function(i) wle.weights(Hmat[,i], y=NULL, smooth=0.0031, 1, raf=1, 
                                                      location=FALSE, max.iter=1000, tol=10^(-6))$weights)
    wH_test = sapply(1:ncol(Hmat_test), function(i) wle.weights(Hmat_test[,i], y=NULL, smooth=0.0031, 1, raf=1, 
                                                                location=FALSE, max.iter=1000, tol=10^(-6))$weights)
    wH = wH/max(wH)
    wH_test = wH_test/max(wH_test)
    wH[wH >= 1-minw] = 1
    wH_test[wH_test >= 1-minw] = 1
    wH[wH <= minw] = 0
    wH_test[wH_test <= minw] = 0
    wH = rowMedians(wH)
    wH_test = rowMedians(wH_test)
    Hmat = sapply(1:ncol(Hmat), function(i) sqrt(wH)*Hmat[,i])
    Hmat_test = sapply(1:ncol(Hmat_test), function(i) sqrt(wH_test)*Hmat_test[,i])
  }
  
  V = NULL
  T = NULL
  T_test = NULL
  
  nbf = dim(Hmat)[2]
  zval = abs(qnorm(alpha/2))
  for(j in 1:hmax){
    try({
      delta = numeric()
      
      for(i in 1:nbf){
        if(model == "ls"){
          log_mod = glm(Y~Hmat[,i], family = binomial)
          logCf = log_mod$coefficients[2]
          sdCf = summary(log_mod)$coefficients[2,2]
          Wval = logCf / sdCf
          if(abs(Wval) <= abs(zval)) logCf = 0
        }else if(model == "wle"){
          log_mod = wle.glm(Y~Hmat[,i], family = binomial)
          logCf = log_mod$root1$coefficients[2]
          sdCf = summary(log_mod)$coefficients[2,2]
          Wval = logCf / sdCf
          if(abs(Wval) <= abs(zval)) logCf = 0
        }
        
        delta[i] = logCf
      }
      delta = as.matrix(delta)
      delta = delta / norm(delta, type = "F")
      
      V = cbind(V, delta)
      T = cbind(T, Hmat %*% delta)
      T_test = cbind(T_test, Hmat_test %*% delta)
      
      for(jj in 1:nbf){
        Hmat[,jj] = Hmat[,jj] - mean(Hmat[,jj]) - as.numeric(cov(Hmat[,jj], T[,j]) / var(T[,j])) * (T[,j] - mean(T[,j]))
        Hmat_test[,jj] = Hmat_test[,jj] - mean(Hmat_test[,jj]) - as.numeric(cov(Hmat_test[,jj], T_test[,j]) / var(T_test[,j])) * (T_test[,j] - mean(T_test[,j]))
      }
      
    }, silent = TRUE)
  }
  
  Vmat = V
  Vmat = Vmat[,!apply(Vmat,2,mean) %in% NaN]
  
  Tmat = T
  Tmat = Tmat[,!apply(Tmat,2,mean) %in% NaN]
  
  Tmat_test = T_test
  Tmat_test = Tmat_test[,!apply(Tmat_test,2,mean) %in% NaN]
  
  if(model == "ls"){
    Fin_log_mod = glm(Y~Tmat, family = binomial)
    Gam_par = Fin_log_mod$coefficients[-1]
    fits = round(Fin_log_mod$fitted.values)
    Gam_constant = Fin_log_mod$coefficients[1]
    
  }else if(model == "wle"){
    Fin_log_mod = wle.glm(Y~Tmat, family = binomial)
    Gam_par = Fin_log_mod$root1$coefficients[-1]
    fits = round(Fin_log_mod$root1$fitted.values)
    Gam_constant = Fin_log_mod$root1$coefficients[1]
  }
  
  preds1 = cbind(1,Tmat_test) %*% as.matrix(c(Gam_constant, Gam_par))
  preds = round(exp(preds1)/(1+exp(preds1)))
  
  Gam_par <- as.matrix(Gam_par)
  
  beta_hat = Vmat%*%Gam_par
  beta_hat_s = Bsf%*%beta_hat
  
  return(list(fits = fits, beta_hat_s = beta_hat_s, beta_hat=beta_hat,
              Gam_constant=Gam_constant, nbf=nbf, preds=preds))
  
}
##____________________________________________________________________________________________________________________________


#_______________Correct classification rate_____________
ccr = function(Y, Ypred){
  n = length(Y)
  p = 0
  for(i in 1:n){
    if(Y[i] == Ypred[i])
      p = p+1
  }
  ccr = p / n
  ce = 1 - ccr
  return(list(ccr = ccr, ce = ce))
}
#_______________________________________________________







