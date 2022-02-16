#' HOMTESTGLM
#'
#' This function tests homogeneity of the regression coefficients in heterogeneous generalized linear models with interactive effects.
#' @param X The (NT) times p design matrix, without an intercept where N=number of individuals, T=length of time series, p=number of explanatory variables.
#' @param Y The T times N panel of response where N=number of individuals, T=length of time series.
#' @param FAMILY A description of the error distribution and link function to be used in the model just like in glm functions.
#' @param Nfactors A pre-specified number of common factors.
#' @param Maxit A maximum number of iterations in optimization. Default is 100.
#' @param tol Tolerance level of convergence. Default is 0.001.
#' @return A list with the following components:
#' \itemize{
#'   \item Coefficients: The estimated heterogeneous coefficients.
#'   \item Factors: The estimated common factors across groups.
#'   \item Loadings: The estimated factor loadings for the common factors.
#'   \item pvalue: The p-value of the homogeneity test.
#' }
#' @references Ando, T. and Bai, J. (2015) A simple new test for slope homogeneity in panel data models with interactive effects. Economics Letters, 136, 112-117.
#' @importFrom stats lm pnorm
#' @export
#' @examples
#' fit <- HOMTESTGLM(data2X,data2Y,binomial(link=logit),2,10,0.5)
HOMTESTGLM<-function(X,Y,FAMILY,Nfactors,Maxit=100,tol=0.001){
  
  
  AY <- Y
  AX <- X
  N <- nrow(Y)
  P <- ncol(Y)
  p <- ncol(X)
  
  
  fit <- PDMIFGLM(AX,AY,binomial(link=logit),2,Maxit,tol)
  B <- (fit$Coefficients)[-1,]
  L <- fit$Loadings
  Fac <- fit$Factors
  Pred <- fit$Predict
  
  H0B <- B%*%rep(1,len=P)/P
  
  ER <- AY-Pred
  
  VAR <- ER^2
  
  LL <- Omega <- matrix(0,p*P,p*P)
  
  A <- matrix(0,p*P,N*Nfactors)
  
  Eta <- matrix(0,p*P,1)
  
  BBB <- matrix(0,p*P,1)
  
  PSI <- matrix(0,Nfactors*N,Nfactors*P)
  
  SSS <- matrix(0,p*P,p*P)
  
  for(j in 1:P){BBB[((j-1)*p+1):(j*p),1] <- B[,j]-H0B}
  
  for(i in 1:N){for(j in 1:P){
    X <- AX[(N*(j-1)+1):(N*j),]; 
    xit <- as.matrix(X[i,],ncol=1)
    lam <- as.matrix(L[j,],ncol=1)
    v <- VAR[i,j]
    a <- v*xit%*%t(lam)
    A[((j-1)*p+1):(j*p),((i-1)*Nfactors+1):(i*Nfactors)] <- a
  }}
  
  ####L
  
  for(i in 1:N){
    for(j in 1:P){
      v <- diag(VAR[i,])
      Psi <- t(L)%*%v%*%L/P
      PSI[((i-1)*Nfactors+1):(i*Nfactors),((j-1)*Nfactors+1):(j*Nfactors)] <- Psi
    }}
  
  
  J <- matrix(0,Nfactors,Nfactors*N)
  for(i in 1:Nfactors){
    v <- as.vector(diag(1,Nfactors)[i,]) 
    v <- rep(v,len= Nfactors*N)
    J[i,] <- v
  }
  
  II <- diagonals::fatdiag(Nfactors*N, steps=N)
  
  for(k in 1:P){
    j <- k
    Ak <- A[((k-1)*p+1):(k*p),]
    Aj <- A[((j-1)*p+1):(j*p),]
    Psi <- PSI[,((j-1)*Nfactors+1):(j*Nfactors)]
    H <- (Psi%*%J)*II
    LL[((j-1)*p+1):(j*p),((k-1)*p+1):(k*p)] <- Ak%*%solve(H)%*%t(Aj)/N
  }
  
  ####SSS
  
  for(j in 1:P){
    v <- diag(VAR[,j])
    X <- AX[(N*(j-1)+1):(N*j),]; 
    SSS[((j-1)*p+1):(j*p),((j-1)*p+1):(j*p)] <- t(X)%*%v%*%X/N
  }
  
  
  ####Eta
  
  for(j in 1:P){
    X <- AX[(N*(j-1)+1):(N*j),]; 
    Eta[((j-1)*p+1):(j*p),1] <- t(X)%*%ER[,j]/sqrt(N)
  }
  
  
  ####Omega
  
  for(j in 1:P){
    X  <- AX[(N*(j-1)+1):(N*j),]; 
    v <- diag(VAR[,j])
    Omega[((j-1)*p+1):(j*p),((j-1)*p+1):(j*p)] <- t(X)%*%v%*%X/N
  }
  
  Swamy <- N*t(BBB)%*%(SSS-LL/P)%*%solve(Omega)%*%t(SSS-LL/P)%*%BBB
  
  z <- as.numeric( (Swamy-P*p)/sqrt(2*P*p) )
  
  pval=2*pnorm(-abs(z))
  
  message("Call:
HOMTESTGLM(X, Y, FAMILY =",FAMILY$family,FAMILY$link,"Nfactors =",Nfactors,", Maxit =",Maxit,", tol =",tol,")
  
N =",P,", T =",N,", p =",p,"

p-value =",pval,"

Fit includes coefficients, factors and loadings.
")
  
  invisible(list("Coefficients"=B,"Factors" = Fac, "Loadings" = L, "pvalue"=pval))
}


