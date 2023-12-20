score.test <- function(model.H0, model.H1) {

 y <- as.matrix(model.frame(model.H0)[,1])
 X <- model.matrix(model.H1)

 coef.post <- which( !(names(coef(model.H1)) %in% names(coef(model.H0))) )
 beta.H0 <- as.matrix(c(coef(model.H1)))
 beta.H0[coef.post] <- 0
 res.H0  <- y-X%*%beta.H0
 sig2.H0 <- sum(res.H0^2)/(length(y)-length(beta.H0))

 U.H0 <- t(X) %*% as.matrix(res.H0)
 I.H0 <- t(X) %*% X
 stat <- t(U.H0) %*% solve(I.H0) %*% U.H0/sig2.H0
 pval <- pchisq(stat,df=1,lower.tail=F)

 out  <- list()
 out[['stat']] <- stat
 out[['pval']] <- pval
out
}


robWald.test <- function(model.H1) {

 y <- as.matrix(model.frame(model.H1)[,1])
 X <- model.matrix(model.H1)
 beta.H1 <- coef(model.H1)
 res.H1  <- as.vector(resid(model.H1))
 sig2.H1 <- sum(res.H1^2)/(length(y)-ncol(X)-1)

 U.H1  <- (X * res.H1)
 rob.vcov <- solve(t(X)%*%X) %*% t(U.H1) %*% U.H1 %*% solve(t(X)%*%X)

 stat <- beta.H1/sqrt(diag(rob.vcov))
 pval <- 2*pnorm(abs(stat),lower.tail=F)
 out  <- list()
 out[['stat']] <- stat
 out[['pval']] <- pval
out
}


LR.test <- function(model.H0,model.H1) {
  stat <- 2*( c(logLik(model.H1)) - c(logLik(model.H0)) )
  pval <- pchisq(stat,lower.tail=F,df=length(coef(model.H1))-length(coef(model.H0)) )
 out  <- list()
 out[['stat']] <- stat
 out[['pval']] <- pval
out
}

