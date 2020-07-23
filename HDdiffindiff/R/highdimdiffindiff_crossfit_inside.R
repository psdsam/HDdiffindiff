#' The core function of doubly robust diff-in-diff estimator under high dimesional covariates
#' @param y0sample1 outcome at time 0 from sample 1
#' @param y1sample1 outcome at time 1 from sample 1
#' @param treatsample1 treatment status at time 1 from sample 1
#' @param xsample1 pretreatment high dimensional covariates from sample 1
#' @param zsample1 pretreatment low diemsioanl covariates allowing nonparametric specification from sample 1
#' @param y0sample2 outcome at time 0 from sample 2
#' @param y1sample2 outcome at time 1 from sample 2
#' @param treatsample2 treatment status at time 1 from sample 2
#' @param xsample2 pretreatment high dimensional covariates from sample 2
#' @param zsample2 pretreatment low diemsioanl covariates allowing nonparametric specification from sample 2
#' @param method basis function in nonparametric specification:
#' \itemize{
#' \item (default) Pol: Polynomial basis with parameter q/2 as the degree of polynomial;
#' \item Tri: Trigonometric polynomials with parameter q/2 as the degree of polynomial;
#' }
#' @param q paramter for the basis function in nonparametric specification
#' @return
#' \item{xdebias}{debiased estimator for the high dimensional coefficients}
#' \item{gdebias}{debiased estimator for the nonparametric coefficients}
#' \item{stdx}{standard error for the high dimensional coefficients}
#' \item{stdg}{standard error for the nonparametric coefficients}
#' @author Sida Peng
#' @export

highdimdiffindiff_crossfit_inside <- function(y0sample1, y1sample1, treatsample1, xsample1, zsample1, y0sample2, y1sample2, treatsample2, xsample2, zsample2, method, q){
  kk1 = dim(xsample1)
  kk2 = dim(xsample2)
  n1 = kk1[1]
  n2 = kk2[1]
  n =  n1+n2
  p = kk1[2]

  pf = c(rep(1, 1,p), rep(0, 1, q))
  if (method=="Pol"){
  zbasissample1 = sieve.Pol(zsample1, q)
  zbasissample2 = sieve.Pol(zsample2, q)
  }else if(method=="Tri"){
    zbasissample1 = sieve.TriPol(zsample1, q/2)
    zbasissample2 = sieve.TriPol(zsample2, q/2)
  }

  wsample1 = cbind(xsample1,zbasissample1[,2:(q+1)])
  wsample2 = cbind(xsample2,zbasissample2[,2:(q+1)])

  fit1 <- cv.glmnet(wsample1,treatsample1,family="binomial",alpha=1,penalty.factor = pf, nfolds = 5)
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
  Phi1fit1 = cv.glmnet(x1,deltay1,alpha=1,penalty.factor = pf)
  Phi0fit1 = cv.glmnet(x0,deltay0,alpha=1,penalty.factor = pf)
  hatPhi1 = predict(Phi1fit1,newx=wsample2,s="lambda.min")
  hatPhi0 = predict(Phi0fit1,newx=wsample2,s="lambda.min")

  newy = rho*(y1sample2-y0sample2-(1-Pi)*hatPhi1 - Pi*hatPhi0)

  newlambda <- cv.glmnet(wsample2,newy,alpha=1, penalty.factor = pf)$lambda.min
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
  biasf = (Sigfinv) %*%Sigf%*%adjy
  gdebias = betahat[(p+1):(p+q)] - biasf

  epsilon = newy - wsample2 %*%betahat #%*%c(xdebias, gdebias)
  Vx = t(covinv)%*% (t(matrix(rep(epsilon, p), ncol = p)*tildex)%*% (matrix(rep(epsilon, p), ncol = p)*tildex)/n) %*%covinv
  Sx = sqrt(diag(Vx))/sqrt(n)
  Of = (t(matrix(rep(epsilon, q), ncol = q)*(zbasissample2[,2:(q+1)]))%*%(matrix(rep(epsilon, q), ncol = q)*(zbasissample2[,2:(q+1)])) -MM%*%t(matrix(rep(epsilon, p), ncol = p)*(xsample2))%*%(matrix(rep(epsilon, p), ncol = p)*(xsample2))%*%t(MM) )/n
  Vf = (Sigfinv*n) %*% Of %*%t(Sigfinv*n)
  Sf = sqrt(diag(Vf))/sqrt(n)
  if (sum(is.na(Sf))>0){
    Vf = (Sigfinv*n) %*%(t(matrix(rep(epsilon, q), ncol = q)*t(Sigf))%*%(matrix(rep(epsilon, q), ncol = q)*t(Sigf))/n) %*%t((Sigfinv*n))
    Sf = sqrt(diag(Vf))/sqrt(n)
  }

  output = list("xdebias" = xdebias,
                "gdebias" = as.vector(gdebias),
                "stdx" = Sx,
                "stdg" = Sf)

  return(output)
}
