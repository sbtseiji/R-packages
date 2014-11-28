cov_select<-function(model=x, cov=y, n.obs=0, AIC=0, cov_orig=NULL){
  require(gym)
  if(is.null(cov_orig)){cov_orig<-cov}
  print(sprintf("AIC= %.4f", AIC ),quote=F)
  #' 偏相関行列を作成
  pmat<- (-1)*solve(cov) / sqrt(diag(solve(cov)) %*% t(diag(solve(cov))))
  diag(pmat)<- 1

  #' 偏相関係数を絶対値に変換
  amat<-abs(pmat)

  #' モデルの係数が0の箇所を無限大に設定
  amat[which(model==0)]<-Inf

  #' 偏相関の絶対値が最小の要素を0にした修正モデルを作成
  model_post<-model
  model_post[which.min(amat)]<-0

  #' モデルのフィットとAICの算出
  f<-fitConGraph(model_post,cov_orig,n.obs)
  AIC_post<-f$dev-2*f$df

  #' モデルの適合度が最大になるまで反復
  if (AIC_post<AIC){
    Recall(model_post,f$Shat,n.obs,AIC=AIC_post,cov_orig=cov)
  }
  
  #' 最終的に得られたモデルを描画 & 偏相関行列を表示
  else{
    drawGraph(model,adjust=FALSE)
    drawGraph(model_post,adjust=FALSE)
    diag(pmat)<-1
 	cat("\n")
    (pmat)
  }
}
