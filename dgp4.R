#__________________data generation process with outliers (Monte Carlo)_______________
library(fda.usc) # Statistical Analysis of Functional Data

dgp4 = function(n1, n2){
  t = seq(1, 21, length.out = 109)
  h = numeric()
  for(i in 1:109)
    h[i] = max(6 - abs(t[i]-11), 0)
  
  h1 = h[51:105]
  h2 = h[1:101]
  h3 = h[9:109]
  
  x1 = matrix(NA, nrow = n1, ncol = 101)
  for(j in 1:n1){
    u = runif(1)
    epsi = rnorm(101, 5, 1)
    x1[j,] = u * h1 + (1-u) * h2 + epsi
  }
  
  x2 = matrix(NA, nrow = n2, ncol = 101)
  for(j in 1:n2){
    u = runif(1)
    epsi = rnorm(101, 5, 1)
    x2[j,] = u * h1 + (1-u) * h3 + epsi
  }
  X = rbind(x1,x2)
  
  return(list(X=X, x1=x1, x2=x2))
}
