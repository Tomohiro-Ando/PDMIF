#' PDMIFLING
#'
#' Under a known group membership, this function estimates heterogeneous panel data models with interactive effects.
#' Together with the regression coefficients, this function estimates the unobserved common factor structures both for across/within groups. 
#' @param X The (NT) times p design matrix, without an intercept where N=number of individuals, T=length of time series, p=number of explanatory variables.
#' @param Y The T times N panel of response where N=number of individuals, T=length of time series.
#' @param Membership A pre-specified group membership.
#' @param NGfactors A pre-specified number of common factors across groups (see example).
#' @param NLfactors A pre-specified number of factors in each groups (see example).
#' @param Maxit A maximum number of iterations in optimization. Default is 100.
#' @param tol Tolerance level of convergence. Default is 0.001.
#' @return A list with the following components:
#' \itemize{
#'   \item Coefficients: The estimated heterogeneous coefficients.
#'   \item Lower05: Lower end (5%) of the 90% confidence interval of the regression coefficients.
#'   \item Upper95: Upper end (95%) of the 90% confidence interval of the regression coefficients.
#'   \item GlobalFactors: The estimated common factors across groups.
#'   \item GlobalLoadings: The estimated factor loadings for the common factors.
#'   \item GroupFactors: The estimated group-specific factors.
#'   \item GroupLoadings: The estimated factor loadings for each group.
#'   \item pval: p-value for testing hypothesis on heterogeneous coefficients.
#'   \item Se: Standard error of the estimated regression coefficients.
#' }
#' @references Ando, T. and Bai, J. (2015) Asset Pricing with a General Multifactor Structure Journal of Financial Econometrics, 13, 556-604. 
#' @importFrom stats lm qnorm pnorm
#' @export
#' @examples
#' N <- 200
#' NGroup <- 3
#' Rs <- rep(2,len=NGroup)
#' PP <- rep(200,len=NGroup)
#' P <- sum(PP)
#' p <- 3
#' R <- 2
#' 
#' AY <- matrix(0,nrow=N,ncol=P)
#' LAB <- sort(rep(1:NGroup,len=P))
#' LAM <- matrix(rnorm(P*R,0,1),nrow=P,ncol=R)
#' FAC <- matrix(rnorm(N*R,0,1),nrow=N,ncol=R)
#' XG <- FAC%*%t(LAM)
#' XL <- matrix(0,ncol=P,nrow=N)
#' 
#' for(i in 1:NGroup){
#'   LAM <- matrix(rnorm(PP[i]*Rs[i],1,1),nrow=PP[i],ncol=Rs[i])
#'   FAC <- matrix(rnorm(N*Rs[i],0,1),nrow=N,ncol=Rs[i])
#'   X <- FAC%*%t(LAM)
#'   XL[,(PP[i]*(i-1)+1):(PP[i]*i)] <- X
#' }
#' 
#' ERR <- matrix(rnorm(N*P,0,1),nrow=N,ncol=P)
#' AY <- XG+XL+ERR
#' AX <- matrix(runif(p*P*N,-2,2),nrow=P*N)
#' AB <- matrix(rnorm(p*P,4,2),ncol=P)
#' for(j in 1:P){AY[,j] <- AY[,j]+AX[(N*(j-1)+1):(N*j),]%*%AB[,j]}
#' 
#' PDMIFLING(AX,AY,LAB,R,Rs)
PDMIFLING <- function (X, Y, Membership, NGfactors, NLfactors, Maxit=100, tol=0.001) 
{
  Ngroups <- length(NLfactors)
  AY <- Y
  AX <- X
  N <- nrow(AY)
  P <- ncol(AY)
  p <- ncol(AX)
  
  PredXB <- matrix(0, nrow = N, ncol = P)
  B <- matrix(0,nrow=(p+1), ncol = P)
  
  for (j in 1:P) {
    X <- AX[(N*(j-1)+1):(N*j),]
    y <- AY[, j]
    fit <- lm(y~X)
    B[, j] <- fit$coefficients
    PredXB[, j] <- cbind(1, X) %*% B[, j]
  }
  
  Z <- AY - PredXB
  VEC <- eigen(Z %*% t(Z))$vectors
  Fac <- sqrt(N) * (VEC)[, 1:NGfactors]
  L <- t(t(Fac) %*% Z/N)
  PredG <- Fac%*%t(L)
  
  Z <- AY - PredXB - PredG
  FS <- matrix(0, nrow = N, ncol = sum(NLfactors))
  LS <- matrix(0, nrow = P, ncol = max(NLfactors))
  PredL <- matrix(0, nrow = N, ncol = P)
  for (i in 1:Ngroups) {
    index <- subset(1:(P), Membership == i)
    Z <- Y[, index]
    VEC <- eigen(Z %*% t(Z))$vectors
    Fac <- sqrt(N) * (VEC)[, 1:NLfactors[i]]
    L <- t(t(Fac) %*% Z/N)
    LS[index, 1:NLfactors[i]] <- L
    if (i == 1) {FS[, 1:NLfactors[1]] <- Fac}
    if (i != 1) {FS[, (sum(NLfactors[1:(i - 1)]) + 1):(sum(NLfactors[1:i]))] <- Fac}
    PredL[, index] <- Fac %*% t(L)
  }
  
  for (ITE in 1:Maxit) {
    B.old <- B
    Y <- AY - PredG - PredL
    for (j in 1:P) {
      X <- AX[(N*(j-1)+1):(N*j),]
      y <- Y[, j]
      fit <- lm(y~X)
      B[, j] <- fit$coefficients
      PredXB[, j] <- cbind(1, X) %*% B[, j]
    }
    
    Z <- AY - PredXB - PredL
    VEC <- eigen(Z %*% t(Z))$vectors
    Fac <- sqrt(N) * (VEC)[, 1:NGfactors]
    L <- t(t(Fac) %*% Z/N)
    PredG <- Fac%*%t(L)
    
    Y <- AY - PredXB - PredG
    for (i in 1:Ngroups) {
      index <- subset(1:(P), Membership == i)
      Z <- Y[, index]
      VEC <- eigen(Z %*% t(Z))$vectors
      Fac <- sqrt(N) * (VEC)[, 1:NLfactors[i]]
      L <- t(t(Fac) %*% Z/N)
      LS[index, 1:NLfactors[i]] <- L
      if (i == 1) {FS[, 1:NLfactors[1]] <- Fac}
      if (i != 1) {FS[, (sum(NLfactors[1:(i - 1)]) + 1):(sum(NLfactors[1:i]))] <- Fac}
      PredL[,index] <- Fac%*%t(L)
    }

    if (mean(abs(B.old - B)) <= tol) {break}
  }
  
  Er <- AY-PredL-PredXB-PredG
  
  Tstat  <- 0*B
  pVal  <- 0*B
  Lower05  <- 0*B
  Upper95  <- 0*B
  V <- 0*B
  B0  <- 0*B #Hypothetical
  S=rep(0,P)
  
  for(i in 1:P){
    S[i] <- mean(Er[,i]^2)
    X <- AX[(N*(i-1)+1):(N*i),]
    W <- cbind(1,X)
    V[,i] <- sqrt(diag(S[i]*t(W)%*%W))
    Lower05[,i]  <- B[,i]+qnorm(0.025)*V[,i]/sqrt(N)
    Upper95[,i]  <- B[,i]+qnorm(0.975)*V[,i]/sqrt(N)
    Tstat[,i] <- sqrt(N)*(B[,i]-B0[,i])/V[,i]
    pVal[,i] <- 2*pnorm(-abs(Tstat[,i]))
  }
  
  return(list("Coefficients"=B,"Lower05"=Lower05,"Upper95"=Upper95,"GlobalFactors"=Fac,
              "GlobalLoadings"=L,"GroupFactors"=FS,"GroupLoadings"=LS,"pval"=pVal,"Se"=V))
  
  
}





