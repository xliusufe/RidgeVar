VAR_MM <- function(y,x,identity=F,Sigma=NULL){
  n = nrow(x)
  p = ncol(x)
  para = c(n,p)

  fit <- .Call("VAR_MM", y,  t(x), Sigma, as.integer(para), as.integer(identity),
               as.integer(is.null(Sigma)))
  fit
}
