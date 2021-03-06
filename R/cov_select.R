cov_select<-function(model=x, cov=y, n.obs=0, AIC=0, cov_orig=NULL){
  require(ggm)
  if(is.null(cov_orig)){cov_orig<-cov}
#  print(sprintf("AIC= %.4f", AIC ),quote=F)
  #' 偏相関行列を作成
  pmat<- (-1)*solve(cov) / sqrt(diag(solve(cov)) %*% t(diag(solve(cov))))
  diag(pmat)<- 1
  
  #' 偏相関係数を絶対値に変換
  amat<-abs(pmat)
  
  #' モデルの係数が0の箇所を無限大に設定
  amat[which(model==0)]<-Inf
  
  #' 偏相関の絶対値が最小の要素を0にした修正モデルを作成
  model_post<-rep(1,nrow(cov))%*%t(rep(1,nrow(cov)))
  model_post[which.min(amat)]<-0
  model_post<-model_post * t(model_post) * model
  
  #' モデルのフィットとAICの算出
  f<-fitConGraph(model_post,cov_orig,n.obs)
  AIC_post<-f$dev-2*f$df
  
  #' モデルの適合度が最大になるまで反復
  if (AIC_post<AIC){
    Recall(model_post,f$Shat,n.obs,AIC=AIC_post,cov_orig=cov_orig)
  }
  
  #' 最終的に得られたモデルを描画 & 偏相関行列を表示
  else{
    diag(pmat)<-1
    pmat[which(model==0)]<-0
    model.0<-model*0
    f<-fitConGraph(model,cov_orig,n.obs)
    f.0<-fitConGraph(model.0,cov_orig,n.obs)

    # GFIの算出
    S<-cov2cor(cov_orig)
    H<-cov2cor(f$Shat)
    p<-nrow(cov_orig)
    num<-sum(diag((solve(H)%*%(S-H))%*%(solve(H)%*%(S-H))))
    den<-sum(diag((solve(H)%*%S)%*%(solve(H)%*%S)))
    GFI<-1-(num/den)
    AGFI<-1-(p*(p+1)*(1-GFI))/(2*f$df)
    RMSEA<-sqrt(max(((f$dev-f$df)/(f$df*(nrow(cov_orig)-1))),0))
    CFI<-(f.0$dev*f.0$df-f$dev*f$df)/(f.0$dev*f.0$df)
    res<-f
    res<-c(f,GFI=GFI, AGFI=AGFI,RMSEA=RMSEA,CFI=CFI)
    return(list(fit=res,model=model, covmat=pmat))
  }
}
