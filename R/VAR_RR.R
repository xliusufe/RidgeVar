VAR_RR <-function(y,x,eta=NULL,alpha=0.1){
  if(is.null(x)) stop("x must not be NA")
  if(is.null(y)) stop("y must not be NA")
  n = nrow(x)
  p = ncol(x)
  q = ncol(y)
  
  if(is.null(p)) p = 1
  if(is.null(q)) q = 1
  is_eta = ifelse(is.null(eta),0,1)
  
  is.debias = F
  para = c(n,p,q,is_eta)
  
  fit <- .Call("VAR_RR", t(x), y, eta, as.integer(para), as.integer(is.debias), alpha)
  fit
}
