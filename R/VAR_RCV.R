VAR_RCV <- function(y,x){
  n = nrow(x)
  p = ncol(x)
  half = ceiling(n/2)
  x1 = x[1:half,]
  y1 = y[1:half]
  x2 = x[-c(1:half),]
  y2 = y[-c(1:half)]

  fit.cv = cv.glmnet(x1,y1,family="gaussian")
  ind = which.min(fit.cv$cvm)
  fits = fit.cv$glmnet.fit
  if(ind==1){
    sigmahat1 = sum(y2^2)/half
  } else{
    ind1 = which(abs(fits$beta[,ind])>0)
    if(length(ind1)>half/2){
      betasort = sort(abs(fits$beta[ind1,ind]),decreasing =TRUE,index.return=T)
      ind1 = ind1[betasort$ix[1:ceiling(half/2)]]
    }
    xm12 = x2[,ind1]
    pm12 = xm12%*%solve(t(xm12)%*%xm12)%*%t(xm12)
    sigmahat1 = (sum(y2^2)-t(y2)%*%pm12%*%y2)/(half-length(ind1))
  }

  fit.cv = cv.glmnet(x2,y2,family="gaussian")
  ind = which.min(fit.cv$cvm)
  fits = fit.cv$glmnet.fit
  if(ind==1){
    sigmahat2 = sum(y1^2)/(n-half)
  } else{
    ind1 = which(abs(fits$beta[,ind])>0)
    if(length(ind1)>half/2){
      betasort = sort(abs(fits$beta[ind1,ind]),decreasing=TRUE,index.return=T)
      ind1 = ind1[betasort$ix[1:ceiling(half/2)]]
    }
    xm21 = x1[,ind1]
    pm21 = xm21%*%solve(t(xm21)%*%xm21)%*%t(xm21)
    sigmahat2 = (sum(y1^2)-t(y1)%*%pm21%*%y1)/(n-half-length(ind1))
  }
  sigmahat = (sigmahat1+sigmahat2)/2
  return(list(sigma2=sigmahat))
}
