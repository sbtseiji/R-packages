partialcor<-function(cormat=x){
  pmat<- (-1)*solve(cormat) / sqrt(diag(solve(cormat)) %*% t(diag(solve(cormat))))
  diag(pmat)<- 1
  (pmat)
}
