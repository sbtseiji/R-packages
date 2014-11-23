drop_var <- function(df = x, min=0.3, max=1,verbose=FALSE,...){
  ncol=length(colnames(df))
  res<-fa(df,...)
  if(verbose){print(res)}
  df2<-df[,apply(abs(res$loadings),1,function(i){
    ifelse((i[order(-i)][1]<min)||(i[order(-i)][2]>=max),
           FALSE,TRUE)})]
  diff=ncol-length(colnames(df2))
  if (diff>0){
    print(paste(diff," item(s) were dropped",sep=""))
    res<-Recall(df2,min=min,max=max,verbose=verbose,...)}
  (res)
}