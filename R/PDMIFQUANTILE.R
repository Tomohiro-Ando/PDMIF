#' PDMIFQUANTILE
#'
#' This function estimates heterogeneous quantile panel data models with interactive effects.
#' @param X The (NT) times p design matrix, without an intercept where N=number of individuals, T=length of time series, p=number of explanatory variables.
#' @param Y The T times N panel of response where N=number of individuals, T=length of time series.
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
#' N <- 300
#' P <- 300
#' R <- 5
#' p <- 8
#' TAU <- 0.95
#' 
#' U <- matrix(runif(N*P,0,1),nrow=N,ncol=P)
#' LAM <- matrix(runif(P*R,-1,1),nrow=P,ncol=R)
#' FAC <- matrix(runif(N*R,0,2),nrow=N,ncol=R)
#' FL <- matrix(0,nrow=N,ncol=P)
#' 
#' for(i in 1:N){
#'   for(j in 1:P){
#'     B <- LAM[j,]+0.1*U[i,j]
#'     if(U[i,j]<=0.2){FL[i,j] <- FAC[i,1:3]%*%B[1:3]}
#'     if(0.2<=U[i,j] && U[i,j]<=0.8){FL[i,j] <- FAC[i,1:4]%*%B[1:4]}
#'     if(0.8<=U[i,j]){FL[i,j] <- FAC[i,1:5]%*%B[1:5]}
#'   }
#' }
#' 
#' AX <- matrix(runif(p*P*N,0,2),nrow=P*N)
#' for(j in 1:P){
#'   AX[(N*(j-1)+1):(N*j),1] <- AX[(N*(j-1)+1):(N*j),1]+0.01*FL[,j]
#'   AX[(N*(j-1)+1):(N*j),3] <- AX[(N*(j-1)+1):(N*j),3]-0.01*FL[,j]
#'   AX[(N*(j-1)+1):(N*j),5] <- AX[(N*(j-1)+1):(N*j),5]+0.02*FL[,j]
#'   AX[(N*(j-1)+1):(N*j),6] <- AX[(N*(j-1)+1):(N*j),6]-0.02*FL[,j]
#' }
#' 
#' XB <- matrix(0,nrow=N,ncol=P)
#' for(i in 1:N){
#'   for(j in 1:P){
#'     X <- AX[(N*(j-1)+1):(N*j),]
#'     B <- c(-1,1,1,-1,rep(1,len=p-4))+0.1*j/N+0.1*U[i,j]
#'     if(U[i,j]<=0.2){XB[i,j] <- X[i,1:p]%*%B[1:p]}
#'     if(0.2<=U[i,j] && U[i,j]<=0.8){XB[i,j] <- X[i,1:p]%*%B[1:p]}
#'     if(0.8<=U[i,j]){XB[i,j] <- X[i,1:p]%*%B[1:p]}
#'   }
#' }
#' 
#' AY <- XB+FL+qnorm(U,0,1)
#' fit <- PDMIFQUANTILE(AX,AY,TAU,R)
PDMIFQUANTILE <- function (X, Y, TAU, Nfactors, Maxit=100, tol=0.001) 
{
  AY <- Y
  AX <- X
  N <- nrow(AY)
  P <- ncol(AY)
  p <- ncol(AX)
  
  
  PredXBB <- matrix(0,nrow=N,ncol=P)
  PredXB <- matrix(0,nrow=N,ncol=P)
  PredFL <- matrix(0,nrow=N,ncol=P)
  
  B <- matrix(0,nrow=p+1,ncol=P)
  
  for(j in 1:P){
    y <- AY[,j]
    X <- AX[(N*(j-1)+1):(N*j),]
    fit <- quantreg::rq(y~X,tau=TAU)
    B[,j] <- (fit$coefficients)
    PredXB[,j] <- cbind(1,X)%*%B[,j]
  }
  
  Z <- AY-PredXB
  VEC <- eigen(Z%*%t(Z))$vectors
  Fac <- sqrt(N)*(VEC)[,1:Nfactors]
  L <- t(solve(t(Fac)%*%Fac)%*%t(Fac)%*%Z)
  PredFL <- Fac%*%t(L)
  
  B.old <- B
  
  for(ITE in 1:Maxit){
    
    for(i in 1:N){
      y <- AY[i,]
      X <- L
      fit <- quantreg::rq(y-PredXB[i,]~X+0,tau=TAU)
      Fac[i,] <- (fit$coefficients)[1:Nfactors]
    }
    
    for(j in 1:P){
      y <- AY[,j]
      X <- cbind(AX[(N*(j-1)+1):(N*j),],Fac)
      fit <- quantreg::rq(y~X,tau=TAU)
      B[,j] <- (fit$coefficients)[1:(p+1)]
      L[j,] <- (fit$coefficients)[(p+2):(p+Nfactors+1)]
      PredXBB[,j] <- AX[(N*(j-1)+1):(N*j),]%*%B[2:(p+1),j]
      PredXB[,j] <- cbind(1,AX[(N*(j-1)+1):(N*j),])%*%B[,j]
      PredFL[,j] <- Fac%*%L[j,]
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
    fit <- summary(quantreg::rq(y~X,tau=TAU),se="iid")
    b <- (fit$coefficients)[1:(p+1),1]
    V[,i] <- (fit$coefficients)[1:(p+1),2]
    Lower05[,i]  <- b+qnorm(0.025)*V[,i]
    Upper95[,i]  <- b+qnorm(0.975)*V[,i]
    pVal[,i] <- (fit$coefficients)[1:(p+1),4]
    Predict[,i] <- quantreg::rq(y~X,tau=TAU)$fitted.values
  }
  
  return(list("Coefficients"=B,"Lower05"=Lower05,"Upper95"=Upper95,
              "Factors"=Fac,"Loadings"=L,"Predict"=Predict,"pval"=pVal,"Se"=V))
}


