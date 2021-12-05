#' PDMIFLIN
#'
#' This function estimates heterogeneous panel data models with interactive effects.
#' This function is similar version of PDMIFLING which accommodates a group structure. 
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
#' @references Ando, T. and Bai, J. (2015) Asset Pricing with a General Multifactor Structure Journal of Financial Econometrics, 13, 556-604. 
#' @importFrom stats lm qnorm
#' @export
#' @examples
#' N <- 100
#' P <- 200
#' 
#' R <- 2
#' p <- 2
#' 
#' LAM <- matrix(rnorm(P*R,0,1),nrow=P,ncol=R)
#' FAC <- matrix(rnorm(N*R,0,1),nrow=N,ncol=R)
#' AY <- FAC%*%t(LAM)+matrix(rnorm(N*P,0,1),nrow=N,ncol=P)
#' 
#' AX <- matrix(0,nrow=N*P,ncol=p)
#' for(j in 1:P){
#'   AX[(N*(j-1)+1):(N*j),1] <- 0.2+0.3*AY[,j]+matrix(rnorm(N,0,1),nrow=N)
#'   AX[(N*(j-1)+1):(N*j),2] <- 0.5+0.5*AY[,j]+matrix(rnorm(N,0,1),nrow=N)
#' }
#' 
#' B <- (1:p)/2-1
#' AB <- B%*%t(rep(1,len=P))
#' 
#' for(j in 1:P){
#' AY[,j] <- AY[,j]+AX[(N*(j-1)+1):(N*j),]%*%AB[,j]
#' }
#' fit <- PDMIFLIN(AX,AY,2)
PDMIFLIN <- function (X, Y, Nfactors, Maxit=100, tol=0.001) 
{
  AX <- X
  AY <- Y
  N <- nrow(AY)
  P <- ncol(AY)
  p <- ncol(AX)
  
  PredXB <- matrix(0, nrow = N, ncol = P)
  B <- matrix(0, nrow=(p+1), ncol = P)
  
  for (j in 1:P) {
    X <- AX[(N*(j-1)+1):(N*j),]
    y <- AY[, j]
    fit <- lm(y~X)
    B[, j] <- fit$coefficients
    PredXB[, j] <- cbind(1, X) %*% B[, j]
  }
  
  Z <- AY - PredXB
  VEC <- eigen(Z %*% t(Z))$vectors
  Fac <- sqrt(N) * (VEC)[, 1:Nfactors]
  L <- t(solve(t(Fac)%*%Fac)%*%t(Fac)%*%Z)
  PredFL <- Fac%*%t(L)
  
  B.old <- B
  
  for (ITE in 1:Maxit) {
    
    Y <- AY - PredFL
    
    for (j in 1:P) {
      X <- AX[(N*(j-1)+1):(N*j),]
      y <- Y[, j]
      fit <- lm(y~X)
      B[, j] <- fit$coefficients
      PredXB[, j] <- cbind(1, X) %*% B[, j]
    }
    
    Z <- AY - PredXB
    VEC <- eigen(Z %*% t(Z))$vectors
    Fac <- sqrt(N) * (VEC)[, 1:Nfactors]
    L <- t(solve(t(Fac)%*%Fac)%*%t(Fac)%*%Z)
    PredFL <- Fac%*%t(L)
    
    if (mean(abs(B.old - B)) <= tol) {
      break
    }
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
    fit <- summary(lm(y~X))
    b <- (fit$coefficients)[1:(p+1),1]
    V[,i] <- (fit$coefficients)[1:(p+1),2]
    Lower05[,i]  <- b+qnorm(0.025)*V[,i]
    Upper95[,i]  <- b+qnorm(0.975)*V[,i]
    pVal[,i] <- (fit$coefficients)[1:(p+1),4]
    Predict[,i] <- lm(y~X)$fitted.values
  }
  
  return(list("Coefficients"=B,"Lower05"=Lower05,"Upper95"=Upper95,
              "Factors" = Fac, "Loadings" = L, "Predict"=Predict, "pval"=pVal,"Se"=V))
  
}
