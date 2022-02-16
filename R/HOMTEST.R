#' HOMTEST
#'
#' This function tests homogeneity of the regression coefficients in heterogeneous panel data models with interactive effects.
#' @param X The (NT) times p design matrix, without an intercept where N=number of individuals, T=length of time series, p=number of explanatory variables.
#' @param Y The T times N panel of response where N=number of individuals, T=length of time series.
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
#' fit <- HOMTEST(data1X,data1Y,2,20,0.5)
HOMTEST<-function(X,Y,Nfactors,Maxit=100,tol=0.001){

  
  AY <- Y
  AX <- X
  N <- nrow(Y)
  P <- ncol(Y)
  p <- ncol(X)
  
  PredXB <- matrix(0,nrow=N,ncol=P)
  B <- matrix(0,nrow=p,ncol=P)
  
  for(j in 1:P){
    y <- AY[,j]
    X <- AX[(N*(j-1)+1):(N*j),]
    fit <- lm(y~X-1)
    B[,j] <- (fit$coefficients)
    PredXB[,j] <- X%*%B[,j]
  }
  
  Z <- AY-PredXB
  VEC <- eigen(Z%*%t(Z))$vectors
  Fac <- sqrt(N)*(VEC)[,1:Nfactors]
  L <- t(solve(t(Fac)%*%Fac)%*%t(Fac)%*%Z)
  PredFL <- Fac%*%t(L)
  
  for(ITE in 1:Maxit){
    
    B.old <- B
    
    Y <- AY-PredFL
    
    for(j in 1:P){
      y <- Y[,j]
      X <- AX[(N*(j-1)+1):(N*j),]
      fit <- lm(y~X-1)
      B[,j] <- (fit$coefficients)
      PredXB[,j] <- X%*%B[,j]
    }
    
    Z <- AY-PredXB
    
    VEC <- eigen(Z%*%t(Z))$vectors
    Fac <- sqrt(N)*(VEC)[,1:Nfactors]
    L <- t(solve(t(Fac)%*%Fac)%*%t(Fac)%*%Z)
    PredFL <- Fac%*%t(L)
    
    if(mean(abs(B.old-B))<=tol){break}
    
  }#ITE
  
  H0B <- B%*%rep(1,len=P)/P
  
  ER <- AY-PredXB-PredFL
  
  S <- W <- Omega <- matrix(0,p*P,p*P)

  M <- diag(1,N)-Fac%*%solve(t(Fac)%*%Fac)%*%t(Fac)
  
  BBB <- matrix(0,p*P,1)
  
  for(j in 1:P){BBB[((j-1)*p+1):(j*p),1] <- B[,j]-H0B}
  
  s <- sum( ER^2 )/(N*P-P*p-(N+P)*Nfactors)
  
  for(j in 1:P){
    
    X <- AX[(N*(j-1)+1):(N*j),]; 
    Omega[((j-1)*p+1):(j*p),((j-1)*p+1):(j*p)] <- s*t(X)%*%M%*%X/N
    S[((j-1)*p+1):(j*p),((j-1)*p+1):(j*p)] <- t(X)%*%M%*%X/N

    for(l in 1:P){
      Xl <- AX[(N*(l-1)+1):(N*l),]; 
      a <- sum( t(L[j,])%*%solve(t(L)%*%L/P)%*%L[l,] )
      W[((j-1)*p+1):(j*p),((l-1)*p+1):(l*p)] <- a*t(X)%*%M%*%Xl/N
    }
    
  }
  
  Swamy <- N*t(BBB)%*%(S-W/P)%*%solve(Omega)%*%t(S-W/P)%*%BBB
  
  z <- as.numeric( (Swamy-P*p)/sqrt(2*P*p) )

  pval=2*pnorm(-abs(z))
  
  message("Call:
HOMTEST(X, Y, Nfactors =",Nfactors,", Maxit =",Maxit,", tol =",tol,")
  
N =",P,", T =",N,", p =",p,"

p-value =",pval,"

Fit includes coefficients, factors and loadings.
")
  
  invisible(list("Coefficients"=B,"Factors" = Fac, "Loadings" = L, "pvalue"=pval))
}


