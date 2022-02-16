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
#' fit <- PDMIFLING(data4X,data4Y,data4LAB,2,c(2,2,2),30,0.1)
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
  Ftemp <- sqrt(N) * (VEC)[, 1:NGfactors]
  Ltemp <- t(t(Ftemp) %*% Z/N)
  FG <- Ftemp
  LG <- Ltemp
  PredG <- FG%*%t(LG)
  
  Z <- AY - PredXB - PredG
  FS <- matrix(0, nrow = N, ncol = sum(NLfactors))
  LS <- matrix(0, nrow = P, ncol = max(NLfactors))
  PredL <- matrix(0, nrow = N, ncol = P)
  for (i in 1:Ngroups) {
    index <- subset(1:(P), Membership == i)
    Z <- Y[, index]
    VEC <- eigen(Z %*% t(Z))$vectors
    Ftemp <- sqrt(N) * (VEC)[, 1:NLfactors[i]]
    Ltemp <- t(t(Ftemp) %*% Z/N)
    LS[index, 1:NLfactors[i]] <- Ltemp
    if (i == 1) {FS[, 1:NLfactors[1]] <- Ftemp}
    if (i != 1) {FS[, (sum(NLfactors[1:(i - 1)]) + 1):(sum(NLfactors[1:i]))] <- Ftemp}
    PredL[, index] <- Ftemp %*% t(Ltemp)
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
    Ftemp <- sqrt(N) * (VEC)[, 1:NGfactors]
    Ltemp <- t(t(Ftemp) %*% Z/N)
    FG <- Ftemp
    LG <- Ltemp
    PredG <- FG%*%t(LG)
    
    Y <- AY - PredXB - PredG
    for (i in 1:Ngroups) {
      index <- subset(1:(P), Membership == i)
      Z <- Y[, index]
      VEC <- eigen(Z %*% t(Z))$vectors
      Ftemp <- sqrt(N) * (VEC)[, 1:NLfactors[i]]
      Ltemp <- t(t(Ftemp) %*% Z/N)
      LS[index, 1:NLfactors[i]] <- Ltemp
      if (i == 1) {FS[, 1:NLfactors[1]] <- Ftemp}
      if (i != 1) {FS[, (sum(NLfactors[1:(i - 1)]) + 1):(sum(NLfactors[1:i]))] <- Ftemp}
      PredL[,index] <- Ftemp%*%t(Ltemp)
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
    Lower05[,i]  <- B[,i]+qnorm(0.05)*V[,i]/sqrt(N)
    Upper95[,i]  <- B[,i]+qnorm(0.95)*V[,i]/sqrt(N)
    Tstat[,i] <- sqrt(N)*(B[,i]-B0[,i])/V[,i]
    pVal[,i] <- 2*pnorm(-abs(Tstat[,i]))
  }
  
  message("Call:
PDMIFLING(X, Y, Memberships, NGfactors =",NGfactors,", NLfactors =",NLfactors,", Maxit =",Maxit,", tol =",tol,")
  
N =",P,", T =",N,", p =",p,"

Fit includes coefficients, confidence interval, global factors, global loadings,
    group factors, group loadings, p-values and standard errors.
")
  
  invisible(list("Coefficients"=B,"Lower05"=Lower05,"Upper95"=Upper95,"GlobalFactors"=FG,
              "GlobalLoadings"=LG,"GroupFactors"=FS,"GroupLoadings"=LS,"pval"=pVal,"Se"=V))
  
  
}





