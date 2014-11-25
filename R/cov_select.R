cov_select<-function(model=x, cov=y, n.obs=0, AIC=Inf){
  require(ggm)

  # 相関行列を作成
  cmat<-cov / sqrt(diag(cov) %*% t(diag(cov)))

  # 偏相関行列を作成
  pmat<- (-1)*solve(cmat) / sqrt(diag(solve(cmat)) %*% t(diag(solve(cmat))))
  diag(pmat)<- 1

  # 相関係数を絶対値に変換
  amat<-abs(cmat)

  # モデルの係数が0の箇所を無限大に設定
  amat[which(model==0)]<-Inf

  # 偏相関の絶対値が最小の要素を0にした修正モデルを作成
  model_post<-model
  model_post[which(amat==min(amat))]<-0

  # モデルのフィットとAICの算出
  f<-fitConGraph(model_post,cov,n.obs)
  AIC_post<-f$dev-2*f$df

  # モデルの適合度が最大になるまで反復
  if (AIC_post<AIC){
    Recall(model_post,f$Shat,n.obs,AIC_post)
  }
  # 最終的に得られたモデルを描画 & 偏相関行列を表示
  else{
    drawGraph(model,adjust=FALSE)
    # 指数表記になるのを防ぐ
    pmat[which(model==0)]<-0
    # 対角成分を1に
    diag(pmat)<-1
    (pmat)
  }
}
