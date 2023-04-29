#__________________data generation process without outliers (Monte Carlo)________________
library(fda.usc) # Statistical Analysis of Functional Data

dgp = function(N){
  dtp = seq(0, 10, length.out = 256)
  Bt = fdata(sin(dtp * pi / 3), dtp)
  nrange = c(0,10)
  n_basis = 13
  
  Bsp = create.bspline.basis(rangeval = nrange, nbasis = n_basis, norder = 4, breaks = 0:10)
  phi_k = eval.basis(dtp, Bsp)
  
  Z = matrix(rnorm(n = (N*n_basis)), ncol = n_basis, byrow = TRUE)
  U = matrix(runif(n = (n_basis*n_basis)), ncol = n_basis, byrow = TRUE)
  
  C_mat = Z %*% U
  Xt = C_mat %*% t(phi_k)
  
  fX=fdata(Xt,dtp)
  temp =  inprod.fdata(fX, Bt)
  
  pihat = exp(temp) / (1 + exp(temp))
  Y = numeric()
  
  for(i in 1:N){
    Y[i] = rbinom(n = 1, size = 1, prob = pihat[i,])
  }
  
  return(list(Y=Y, X=Xt, Bt=Bt))
}

