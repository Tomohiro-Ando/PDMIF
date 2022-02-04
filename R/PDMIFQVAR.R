#' PDMIFQVAR
#'
#' This function estimates heterogeneous quantile panel data VAR models with interactive effects.
#' @param Y The T times N panel of response where N=number of individuals, T=length of time series.
#' @param LAG The number of lags from y_t-1 to y_t-LAG used in the VAR.
#' @param TAU A pre-specified quantile point.
#' @param Nfactors A pre-specified number of common factors.
#' @param Maxit A maximum number of iterations in optimization. Default is 100.
#' @param tol Tolerance level of convergence. Default is 0.001.
#' @return A list with the following components:
#' \itemize{
#'   \item Coefficients: The estimated heterogeneous coefficients.
#'   \item Lower05: Lower end (5%) of the 90% confidence interval of the regression coefficients.
#'   \item Upper95: Upper end (95%) of the 90% confidence interval of the regression coefficients.
#'   \item Factors: The estimated common factors across groups.
#'   \item Loadings: The estimated quantile point under a given tau.
#'   \item Predict: The conditional expectation of response variable. 
#'   \item pval: p-value for testing hypothesis on heterogeneous coefficients.
#'   \item Se: Standard error of the estimated regression coefficients.
#' }
#' @references Ando, T. and Bai, J. (2020) Quantile co-movement in financial markets Journal of the American Statistical Association. 
#' @importFrom stats qnorm
#' @export
#' @examples
#' fit <- PDMIFQVAR(data8Y,2,0.1,2)
PDMIFQVAR <- function (Y, LAG, TAU, Nfactors, Maxit=100, tol=0.001) 
{

  AX <- Y[LAG:(nrow(Y)-1),]
  if (LAG>1){
  for(i in (1:(LAG-1))){
   AX <- cbind(AX,Y[(LAG-i):(nrow(Y)-1-i),]) 
  }  }
  AY <- Y[(LAG+1):(nrow(Y)),]
  
N <- nrow(AY)
P <- ncol(AY)
p <- ncol(AX)

if(p+Nfactors<N){
  
  XB <- matrix(0,nrow=N,ncol=P)
  FL <- matrix(0,nrow=N,ncol=P)
  B <- matrix(0,nrow=p+1,ncol=P)
  
  for(j in 1:P){
    y <- AY[,j]
    X <- AX
    fit <- quantreg::rq(y~X,tau=TAU)
    B[,j] <- (fit$coefficients)
    XB[,j] <- cbind(1,X)%*%B[,j]
  }
  
  Z <- AY-XB
  VEC <- eigen(Z%*%t(Z))$vectors; 
  Fac <- sqrt(N)*(VEC)[,1:Nfactors]; 
  L <- t(solve(t(Fac)%*%Fac)%*%t(Fac)%*%Z);
  FL <- Fac%*%t(L)
  
  B.old <- B
  FL.old <- FL
  
  for(ITE in 1:Maxit){
    
    for(i in 1:N){
      y <- AY[i,]
      X <- L
      fit <- quantreg::rq(y-XB[i,]~X+0,tau=TAU)
      Fac <- as.matrix(Fac)
      Fac[i,] <- (fit$coefficients)[1:Nfactors]
    }
    
    for(j in 1:P){
      y <- AY[,j]
      X <- cbind(AX,Fac)
      fit <- summary(quantreg::rq(y~X,tau=TAU))
      B[,j] <- (fit$coefficients)[1:(p+1),1]
      L[j,] <- (fit$coefficients)[(p+2):(p+Nfactors+1),1]
      XB[,j] <- cbind(1,AX)%*%B[,j]
      FL[,j] <- Fac%*%L[j,]
    }
    
    Up <- sum((B.old-B)^2)/(length(B))+sum((FL.old-FL)^2)/(length(FL))
    
    B.old <- B
    FL.old <- FL
    
    if(Up<=tol){break}
    
  } # close ITE loop
}

  pVal  <- 0*B
  Lower05  <- 0*B
  Upper95  <- 0*B
  V  <- 0*B
  Predict <- 0*AY
  for(i in 1:P){
    y <- AY[,i]
    X <- cbind(AX,Fac)
    fit <- summary(quantreg::rq(y~X,tau=TAU),se="iid")
    b <- (fit$coefficients)[1:(p+1),1]
    V[,i] <- (fit$coefficients)[1:(p+1),2]
    Lower05[,i]  <- b+qnorm(0.05)*V[,i]
    Upper95[,i]  <- b+qnorm(0.95)*V[,i]
    pVal[,i] <- (fit$coefficients)[1:(p+1),4]
    Predict[,i] <- quantreg::rq(y~X,tau=TAU)$fitted.values
  }

cat("Call:
PDMIFQVAR(Y, LAG =",LAG,", TAU =",TAU,", Nfactors =",Nfactors,", Maxit =",Maxit,", tol =",tol,")
  
N =",P,", T =",N,"

Fit includes coefficients, confidence interval, factors, loadings,
    expected values, p-values and standard errors.
")


invisible(list("Coefficients"=B,"Lower05"=Lower05,"Upper95"=Upper95,
               "Factors"=Fac,"Loadings"=L,"Predict"=Predict,"pval"=pVal,"Se"=V))
}
