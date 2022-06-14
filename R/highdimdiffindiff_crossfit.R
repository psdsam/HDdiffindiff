highdimdiffindiff_crossfit3 <- function(y0, y1, treat, x, z, q, k, z0, alp){
  source("highdimdiffindiff_crossfit_inside3.r")
  library(glmnet)
  library(MASS)
  library(flare)
  kk = dim(x)
  n0 = kk[1]
  p = kk[2]
  qq = length(z0)
  idx = c(1:n0)
  group = ceiling(idx/ceiling(n0/k))
  xdebias = rep(0, p)
  gdebias = rep(0, qq)
  stdx = rep(0, p)
  stdg = rep(0, qq)
  for (i in 1:k){
     zsample1 = z[group!=i]
     xsample1 = x[group!=i,]
     y0sample1 = y0[group!=i]
     y1sample1 = y1[group!=i]
     treatsample1 = treat[group!=i]

     zsample2 = z[group==i]
     xsample2 = x[group==i,]
     y0sample2 = y0[group==i]
     y1sample2 = y1[group==i]
     treatsample2 = treat[group==i]
     
     ff = highdimdiffindiff_crossfit_inside3(y0sample1, y1sample1, treatsample1, xsample1, zsample1, y0sample2, y1sample2, treatsample2, xsample2, zsample2, q, z0, alp)
     xdebias = xdebias+ff$xdebias
     gdebias = gdebias + ff$gdebias
     stdx = stdx+ff$stdx^2
     stdg = stdg+ff$stdg^2
  }
  xdebias = xdebias/k
  gdebias = gdebias/k
  stdx = sqrt(stdx/(k^2))
  stdg = sqrt(stdg/(k^2))
  CIpoint = rbind(c(xdebias, gdebias)+qnorm(alp/2)*c(stdx, stdg), c(xdebias, gdebias)+qnorm(1-alp/2)*c(stdx, stdg))
  CIuniform = rbind(debias+ff$tc[1]*stdg, gdebias+ff$tc[2]*stdg)
  #browser()
  output = list("xdebias" = xdebias,
                "gdebias" = gdebias,
                "stdx" =  stdx,
                "stdg" = stdg,
                "CIpoint" = CIpoint,
                "CIuniform" = CIuniform)
  
  return(output)
  
}
  
  
  