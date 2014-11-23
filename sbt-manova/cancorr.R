# 正準相関分析の結果要約
# object　cancorrクラスオブジェクト

summary.cancorr<-function (object)
{
    ans<-list("Canonical Correlations"= object$cor, 
              "Eigenvalues of CanRsq/(1-CanRsq)"=object$values,
              "Test of Canonical Correlations"=object$statistics,
              "Cannonical Coefficients for X variables"=object$xcoeff,
              "Cannonical Coefficients for Y variables"=object$ycoeff)
    ans$"Canonical Structure"<-noquote("and their Canonical Variables   ")
    names(ans$"Canonical Structure")<-"Correlations between Raw Variables"    
    ans$"X variables and their Canonical Variables" <-object$xcor
    ans$"Y variables and their Canonical Variables" <-object$ycor
    if(!is.null(object$xpro)){

        ans$"Canonical Redundancy Analysis"<-
              noquote("Canonical Variables  <--------->   Canonical Variables")
        names(ans$"Canonical Redundancy Analysis")<- 
              "Their Own                          The Opposite       "
        xpsum <-cumsum(object$xpro)
        xrsum <-cumsum(object$xred)
        ans$"X Variables"<-data.frame(
              "Proportion"=object$xpro,
              "Cumulative"=xpsum,
              "Can.R-Sq."=object$cor^2,
              "Proportion"=object$xred,
              "Cumulative"=xrsum)

        ypsum <-cumsum(object$ypro)
        yrsum <-cumsum(object$yred)
        ans$"Y Variables"<-data.frame(
              "Proportion"=object$ypro,
              "Cumulative"=ypsum,
              "Can.R-Sq."=object$cor^2,
              "Proportion"=object$yred,
              "Cumulative"=yrsum)
    }
    
    if(!is.null(object$score)) ans$"Canonical Scores"<-object$score
    
    class(ans) <-c("summary.cancorr","listof")
    ans
}



# 正準相関分析
# x　データ行列またはデータフレーム
# y　データ行列またはデータフレーム
# scale　TRUE のときデータの正規化を行う (論理値)
# test　正準変量の検定方法 (”ChiSquare”または”Wilks”)
# redundancy　TRUE のとき冗長性分析を行う (論理値)
# score　TRUE のとき正準得点の算出を行う (論理値)

cancorr<-function(x,y,scale=TRUE, test="Wilks",redundancy=TRUE,score=FALSE) 
{
    if(scale){
        xx <- scale(as.matrix(x))
        yy <- scale(as.matrix(y))
    }
    else{
        xx <- as.matrix(x)
        yy <- as.matrix(y)
    }
	tmp<- cbind(xx,yy)
	tmp<- subset(tmp,complete.cases(tmp))
	nx<- ncol(xx)
	ny<- ncol(xx) +1
	xx<- tmp[,1:nx]
	yy<- tmp[,ny:ncol(tmp)]
    n  <- nrow(xx)
    cor <-cancor(xx,yy,xcenter=F,ycenter=F)
    xcoef <- sqrt(n-1)*cor$xcoef
    ycoef <- sqrt(n-1)*cor$ycoef
    V <- xx %*% xcoef
    W <- yy %*% ycoef

    VX <- t(xx)%*%V/(n-1)
    WY <- t(yy)%*%W/(n-1)
    U  <- -log(1-cor$cor^2)
    k  <- seq(along=cor$cor)-1
    m <- NROW(cor$cor)
    e <- cor$cor^2/(1-cor$cor^2)
    prop <- e/sum(e)
    cum <- cumsum(prop)
    emat <-data.frame(EigenValue=e, Proportion=prop, Cumulative=cum)
    
    if(test=="ChiSquare"){
        df <- (nrow(xcoef)-k)*(nrow(ycoef)-k)
        Chisq <- (nrow(V)-(nrow(xcoef)+nrow(ycoef)+1)/2)*rev(cumsum(rev(U)))
        pval <- pchisq(Chisq,df,lower=FALSE)
        STAT <-data.frame(ChiSq=Chisq, df=df, Prob=round(pval,8))
    }
    else{
        df<-max(ncol(xx),ncol(yy))
        STAT<-matrix(c(rep(0,NROW(e)*5)),ncol=5)
        for(i in 1:NROW(e)){
            lambda <- Wilks(e[i:NROW(e)],df-i+1,n-df-1)
            q <- lambda[2]
            df1 <- lambda[3]
            df2 <- lambda[4]
            pval <- pf(q,df1,df2,lower=F)
            STAT[i,] <- c(lambda[1],q,df1,df2,round(pval,6))
        }
        colnames(STAT) <-c("Lambda","ApproxF","NumDF","DenDF","Prob")
        STAT<-as.data.frame(STAT)
    }
    
    colnames(xcoef) <-colnames(xcoef,do.NULL=FALSE,prefix="V")
    rownames(xcoef) <-colnames(x,do.NULL=FALSE)
    colnames(ycoef) <-colnames(ycoef,do.NULL=FALSE,prefix="W")
    rownames(ycoef) <-colnames(y,do.NULL=FALSE)
    colnames(VX) <-colnames(VX,do.NULL=FALSE,prefix="V")
    colnames(WY) <-colnames(WY,do.NULL=FALSE,prefix="W")

    xpro<-ypro<-xred<-yred<-cscore<-NULL
    
    if(redundancy){
       xpro  <-apply(VX^2,2,sum)/nrow(VX)
       ypro  <-apply(WY^2,2,sum)/nrow(WY)
       xred  <-xpro[1:NROW(cor$cor)]*(cor$cor^2)
       yred  <-ypro[1:NROW(cor$cor)]*(cor$cor^2)
    }

    if(score){
        xs <-xx%*%xcoef[,1:m]
        ys <-yy%*%ycoef[,1:m]

        cscore <- cbind(xs,ys)
        rownames(cscore) <-rownames(cscore,do.NULL=FALSE,prefix="Obs")
    }

    structure(list("cor" = cor$cor,
                   "values"=emat,
                   "statistics"=STAT,
                   "xcoeff" = xcoef[,1:m],
                   "ycoeff" = ycoef[,1:m],
                   "xcor"=VX[,1:m],
                   "ycor"=WY[,1:m],
                   "xpro"=xpro[1:m],
                   "ypro"=ypro[1:m],
                   "xred"=xred,
                   "yred"=yred,
                   "score"=cscore),
                   class="cancorr")
 }