


highdimdiffindiff_crossfit_inside3 <- function(y0sample1, y1sample1, treatsample1, xsample1, zsample1, y0sample2, y1sample2, treatsample2, xsample2, zsample2, q, z0, alp){
  source("Sieve_Functional_Space_Basis.r")
  kk1 = dim(xsample1)
  kk2 = dim(xsample2)
  n1 = kk1[1]
  n2 = kk2[1]
  n =  n1+n2
  p = kk1[2]
  qq = length(z0)
  # zsample1 = z[1:floor(n0*9/10)]
  # xsample1 = x[1:floor(n0*9/10),]
  # y0sample1 = y0[1:floor(n0*9/10)]
  # y1sample1 = y1[1:floor(n0*9/10)]
  # treatsample1 = treat[1:floor(n0*9/10)]
  
  # zsample1 = z[(floor(n0/2)+1):n0]
  # xsample1 = x[(floor(n0/2)+1):n0,]
  # y0sample1 = y0[(floor(n0/2)+1):n0]
  # y1sample1 = y1[(floor(n0/2)+1):n0]
  # treatsample1 = treat[(floor(n0/2)+1):n0]
  # zsample1 = z[(floor(n0/2)+1):n0]
  
  # xsample2 = x[(floor(n0*9/10)+1):n0,]
  # y0sample2 = y0[(floor(n0*9/10)+1):n0]
  # y1sample2 = y1[(floor(n0*9/10)+1):n0]
  # treatsample2 = treat[(floor(n0*9/10)+1):n0]
  # zsample2 = z[(floor(n0*9/10)+1):n0]
  

  #browser()
  pf = c(rep(1, 1,p), rep(0, 1, q))
  zbasissample1 = sieve.Pol(zsample1, q) #sieve.TriPol(zsample1, q/2) 
  zbasissample2 = sieve.Pol(zsample2, q) #sieve.TriPol(zsample2, q/2)
  
  zbasispredict = sieve.Pol(z0, q) #sieve.TriPol(z0, q/2)
  zbasispredict = zbasispredict[,2:(q+1)]
  
  wsample1 = cbind(xsample1,zbasissample1[,2:(q+1)])
  wsample2 = cbind(xsample2,zbasissample2[,2:(q+1)])

  #fit1 <- cv.glmnet(wsample1,treatsample1,family="binomial",alpha=1,penalty.factor = pf, nfolds = 5)
  fit1 <- cv.glmnet(wsample1,treatsample1,family="binomial",alpha=1,nfolds = 3)
  Pi = predict(fit1,newx=wsample2,type='response',s="lambda.min")

  #delete data with 0 propensity score
  isnanidx = which(Pi>0.99|Pi<0.01)
  if (length(isnanidx)!=0){
    y1sample2 = y1sample2[-isnanidx]
    y0sample2 = y0sample2[-isnanidx]
    treatsample2 = treatsample2[-isnanidx]
    Pi = Pi[-isnanidx]
    wsample2 = wsample2[-isnanidx,]
    xsample2 = xsample2[-isnanidx,]
    zbasissample2 = zbasissample2[-isnanidx,]
  }
  kk = dim(xsample2)
  n = kk[1]
  p = kk[2]
  
  rho = (treatsample2-Pi)/(Pi*(1-Pi))
  
  x1 = wsample1[treatsample1==1,]
  x0 = wsample1[treatsample1==0,]
  deltay1 = y1sample1[treatsample1==1]-y0sample1[treatsample1==1]
  deltay0 = y1sample1[treatsample1==0]-y0sample1[treatsample1==0]
  #Phi1fit1 = cv.glmnet(x1,deltay1,alpha=1,penalty.factor = pf)
  #Phi0fit1 = cv.glmnet(x0,deltay0,alpha=1,penalty.factor = pf)
  Phi1fit1 = cv.glmnet(x1,deltay1,alpha=1, nfold=5)
  Phi0fit1 = cv.glmnet(x0,deltay0,alpha=1, nfold=5)
  hatPhi1 = predict(Phi1fit1,newx=wsample2,s="lambda.min")
  hatPhi0 = predict(Phi0fit1,newx=wsample2,s="lambda.min")
  
  newy = rho*(y1sample2-y0sample2-(1-Pi)*hatPhi1 - Pi*hatPhi0)
  
  newlambda <- cv.glmnet(wsample2,newy,alpha=1, penalty.factor = pf, nfold=3)$lambda.min
  fit2 <- glmnet(wsample2,newy,alpha=1, penalty.factor = pf, lambda = newlambda)
  betahat <- fit2$beta
  
  #debias
  tildex = xsample2 - zbasissample2%*%solve(t(zbasissample2)%*%zbasissample2)%*%t(zbasissample2)%*%xsample2
  out1 =  sugm(tildex, method = "clime")
  out = sugm.select(out1, criterion ="cv", loss = "tracel2")
  covinv = matrix(unlist(out$opt.icov), ncol = p, byrow = TRUE)
  
  adjy = (rho*(y1sample2-y0sample2-(1-Pi)*hatPhi1 - Pi*hatPhi0 )- predict(fit2, newx = wsample2))
  xdebias = betahat[1:p] + t(adjy)%*%tildex%*%covinv/n
  
  hatM = matrix(0, q, n)
  MM = matrix(0, q,p)
  for (ii in 1:q){
    out2 =  cv.glmnet(xsample2, zbasissample2[,ii+1],alpha=1, intercept=FALSE)
    MM[ii,] = coef(out2)[2:(p+1)]
    hatM[ii,] = t(zbasissample2[,ii+1] - predict(out2,newx=xsample2,s="lambda.min"))
  }
  
  Sigf = t(zbasissample2[,2:(q+1)] - xsample2%*%t(MM))
  Sigfinv = solve(Sigf%*%zbasissample2[,2:(q+1)])
  biasf = (Sigfinv) %*%Sigf%*%adjy #(Sigfinv%*%hatM)%*%adjy
  gdebias = (zbasispredict)%*%(betahat[(p+1):(p+q)] - biasf)
  browser
  #if (abs(gdebias[1])>2){
  #   browser()
  # }
  
  epsilon = newy - wsample2 %*%betahat #%*%c(xdebias, gdebias)
  Vx = t(covinv)%*% (t(matrix(rep(epsilon, p), ncol = p)*tildex)%*% (matrix(rep(epsilon, p), ncol = p)*tildex)/n) %*%covinv
  Sx = sqrt(diag(Vx))/sqrt(n)
  #browser()
  Of = (t(matrix(rep(epsilon, q), ncol = q)*(zbasissample2[,2:(q+1)]))%*%(matrix(rep(epsilon, q), ncol = q)*(zbasissample2[,2:(q+1)])) -MM%*%t(matrix(rep(epsilon, p), ncol = p)*(xsample2))%*%(matrix(rep(epsilon, p), ncol = p)*(xsample2))%*%t(MM) )/n
  Vf = (Sigfinv*n) %*% Of %*%t(Sigfinv*n)
  Sf = sqrt(diag((zbasispredict)%*%Vf%*%t(zbasispredict)))/sqrt(n)
  if (sum(is.na(Sf))>0){
    Vf = (Sigfinv*n) %*%(t(matrix(rep(epsilon, q), ncol = q)*t(Sigf))%*%(matrix(rep(epsilon, q), ncol = q)*t(Sigf))/n) %*%t((Sigfinv*n))
    Sf = sqrt((zbasispredict)%*%Vf%*%t(zbasispredict))/sqrt(n)
    #Sf = sqrt(diag(Vf))/sqrt(n)
  }
  sfs = matrix(Sf, nrow=1000, ncol=qq, byrow=TRUE)
  Vfdecp = svd(Vf)
  bootrand = matrix(rnorm(1000*q,0, 1),nrow=q)
  tboot = (zbasispredict)%*%(Vfdecp$u%*%diag(Vfdecp$d^(1/2))%*%t(Vfdecp$v))%*%bootrand/sqrt(n)/t(sfs)     
  tc =  apply(tboot, 1, quantile, probs = c(alp/2, 1-alp/2),  na.rm = TRUE)
  tc = c(min(tc[1,]), max(tc[2,]))
  #browser()
  output = list("xdebias" = xdebias,
                "gdebias" = as.vector(gdebias),
                "stdx" = Sx,
                "stdg" = Sf,
                "tc" = tc)
  
  return(output)
} 