#__________________data generation process with outliers (Monte Carlo)________________
library(fda.usc) # Statistical Analysis of Functional Data

  dgp2 = function(N){
  dtp = seq(0, 10, length.out = 256)
  Bt = fdata(2 * sin(dtp * 4*pi / 3), dtp)
  
  nrange = c(0,10)
  n_basis = 13
  
  Bsp = create.bspline.basis(rangeval = nrange, nbasis = n_basis, norder = 4, breaks = 0:10)
  phi_k = eval.basis(dtp, Bsp)
  
  Z = matrix(rnorm(n = (N*n_basis)), ncol = n_basis, byrow = TRUE)
  U = matrix(runif(n = (n_basis*n_basis), min = 0, max = 5), ncol = n_basis, byrow = TRUE)
  
  C_mat = Z %*% U
  Xt = C_mat %*% t(phi_k)
  
  fX=fdata(Xt,dtp)
  temp =  inprod.fdata(fX, Bt)
  
  pihat = exp(temp) / (1 + exp(temp))
  Y = numeric()
  
  for(i in 1:N){
    Y[i] = rbinom(n = 1, size = 1, prob = pihat[i,])
  }
  return(list(X = Xt, Y = Y, Bt=Bt))
}



