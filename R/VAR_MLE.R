VAR_MLE <- function(y,x,max.iter=50,tol=1e-4){
  n = nrow(x)
  p = ncol(x)
  para = c(n,p)

  fit <- .Call("VAR_MLE", y, t(x), as.integer(para), as.integer(max.iter),
               as.double(tol))
  fit
}
