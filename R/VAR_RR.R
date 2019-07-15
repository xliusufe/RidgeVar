VAR_RR <-function(y,x,eta=NULL,alpha=0.1){
  n = nrow(x)
  p = ncol(x)
  para = c(n,p)
  is.debias=F

  fit <- .Call("VAR_RR", y, t(x), eta, as.integer(para), as.integer(is.debias),
               alpha, as.integer(is.null(eta)))
  fit
}
