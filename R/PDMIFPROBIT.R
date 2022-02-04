#' PDMIFPROBIT
#'
#' This function estimates heterogeneous probit panel data models with interactive effects.
#' @param X The (NT) times p design matrix, without an intercept where N=number of individuals, T=length of time series, p=number of explanatory variables.
#' @param Y The T times N panel of response where N=number of individuals, T=length of time series.
#' @param Nfactors A pre-specified number of common factors.
#' @param Maxit A maximum number of iterations in optimization. Default is 100.
#' @param tol Tolerance level of convergence. Default is 0.001.
#' @return A list with the following components:
#' \itemize{
#'   \item Coefficients: The estimated heterogeneous coefficients.
#'   \item Lower05: Lower end (5%) of the 90% confidence interval of the regression coefficients.
#'   \item Upper95: Upper end (95%) of the 90% confidence interval of the regression coefficients.
#'   \item Factors: The estimated common factors across groups.
#'   \item Loadings: The estimated factor loadings for the common factors.
#'   \item Predict: The conditional expectation of response variable. 
#'   \item pval: p-value for testing hypothesis on heterogeneous coefficients.
#'   \item Se: Standard error of the estimated regression coefficients.
#' }
#' @references Ando, T., Bai, J. and Li, K. (2021) Bayesian and maximum likelihood analysis of large-scale panel choice models with unobserved heterogeneity, Journal of Econometrics. 
#' @importFrom stats glm qnorm binomial
#' @export
#' @examples 
#' fit <- PDMIFPROBIT(data2X,data2Y,2)
PDMIFPROBIT <- function (X, Y, Nfactors, Maxit=100, tol=0.001) 
{
  AY <- Y
  AX <- X
  N <- nrow(AY)
  P <- ncol(AY)
  p <- ncol(AX)
  
  PredXB <- matrix(0,nrow=N,ncol=P)
  B <- matrix(0,nrow=p+1,ncol=P)
  
  for(j in 1:P){
    y <- AY[,j]
    X <- AX[(N*(j-1)+1):(N*j),]
    fit <- glm(y~X,family=binomial(link=probit))
    B[,j] <- (fit$coefficients)
    PredXB[,j] <- cbind(1,X)%*%B[,j]
  }
  
  Z <- AY-PredXB
  VEC <- eigen(Z%*%t(Z))$vectors
  Fac <- (VEC)[,1:Nfactors]
  L <- t(solve(t(Fac)%*%Fac)%*%t(Fac)%*%Z)
  PredFL <- Fac%*%t(L)
  
  AB <- rbind(B,t(L))
  B.old <- B
  
  for(ITE in 1:Maxit){
    
    for(j in 1:P){
      y <- AY[,j]
      X <- cbind(AX[(N*(j-1)+1):(N*j),],Fac)
      fit <- glm(y~X,family=binomial(link=probit))
      AB[,j] <- fit$coefficients
      L[j,] <- AB[-(1:(p+1)),j]
      B[,j] <- AB[1:(p+1),j]
      PredXB[,j] <- cbind(1,AX[(N*(j-1)+1):(N*j),])%*%AB[1:(p+1),j]
    }
    
    for(j in 1:N){
      y <- AY[j,]
      X <- L
      fit <- glm(y~X+offset(PredXB[j,])+0,family=binomial(link=probit))
      if(Nfactors==1){Fac[j] <- fit$coefficients}
      if(Nfactors!=1){Fac[j,] <- fit$coefficients}
    }
    
    if(mean(abs(B.old-B))<=tol){break}
    
    B.old <- B
    
  }#ITE
  
  pVal  <- 0*B
  Lower05  <- 0*B
  Upper95  <- 0*B
  V  <- 0*B
  Predict <- 0*AY

  for(i in 1:P){
    y <- AY[,i]
    X <- cbind(AX[(N*(i-1)+1):(N*i),],Fac)
    fit <- summary(glm(y~X,family=binomial(link=probit)))
    b <- (fit$coefficients)[1:(p+1),1]
    V[,i] <- (fit$coefficients)[1:(p+1),2]
    Lower05[,i]  <- b+qnorm(0.05)*V[,i]
    Upper95[,i]  <- b+qnorm(0.95)*V[,i]
    pVal[,i] <- (fit$coefficients)[1:(p+1),4]
    Predict[,i] <- glm(y~X,family=binomial(link=probit))$fitted.values
  }
  
  cat("Call:
PDMIFPROBIT(X, Y, Nfactors =",Nfactors,", Maxit =",Maxit,", tol =",tol,")
  
N =",P,", T =",N,", p =",p,"

Fit includes coefficients, confidence interval, factors, loadings,
    expected values, p-values and standard errors.
")
  
  invisible(list("Coefficients"=B,"Lower05"=Lower05,"Upper95"=Upper95,
              "Factors"=Fac,"Loadings"=L,"Predict"=Predict,"pval"=pVal,"Se"=V))
}
