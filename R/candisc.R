candisc<-function(object)  #多変量 lmオブジェクト
{
    require(MASS)
  
    m  <- min(object$rank-1,ncol(object$coefficients))

    w  <- resid(object)
    W  <- t(w)%*%w
    t  <- resid(update(object,~1))

    T  <- t(t)%*%t 
    B  <- T-W
    We <- eigen(W,symmetric=TRUE)
    W1 <- solve(t(We$vectors %*% diag(sqrt(We$values))))
    ev <- eigen( t(W1)%*% B %*% W1)
    cv <- W1 %*% ev$vectors * sqrt(nrow(object$model)-object$rank)
    cv <- cv[,1:m]

    e <- ev$values
    prop <- e[1:m]/sum(e[1:m])
    cum <- cumsum(prop)
    emat <-data.frame(EigenValue=e[1:m], Proportion=prop, Cumulative=cum)
   
    if (m>1){
        colnames(cv) <-colnames(cv,do.NULL=FALSE,prefix="CAN")
        rownames(cv) <-colnames(object$coefficients,do.NULL=FALSE)
    }
    
    STAT<-matrix(c(rep(0,m*5)),ncol=5)
    for(i in 1:m){
        lambda <- Wilks(e[i:NROW(e)],object$rank-i,object$df.residual)
        q      <- lambda[2]
        df1    <- lambda[3]
        df2    <- lambda[4]
        pval   <- round(pf(q,df1,df2,lower=F),8)
        STAT[i,] <- c(lambda[1],q,df1,df2,pval)
    }
    colnames(STAT) <-c("Lambda","ApproxF","NumDF","DenDF","Prob")
    STAT<-as.data.frame(STAT)
    
    structure(list( "EigenValues"=emat,
                    "Statistics"=STAT,
                    "Canonical Variate Coefficients"=cv))    
}