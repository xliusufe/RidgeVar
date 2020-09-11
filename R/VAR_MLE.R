VAR_MLE <- function(y,x,max.iter=50,tol=1e-4){
  if(is.null(x)) stop("x must not be NA")
  if(is.null(y)) stop("y must not be NA")
  n = nrow(x)
  p = ncol(x)
  if(is.null(p)) p = 1
  para = c(n,p)

  fit <- .Call("VAR_MLE", y, t(x), as.integer(para), as.integer(max.iter), as.double(tol))
  fit
}
