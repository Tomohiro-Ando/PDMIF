#' PDMIFCLUST
#'
#' Under a pre-specified number of groups and the number of common factors, this function implements clustering for N individuals in the panels.
#' Each of individuals in the group are subject to the group-specific unobserved common factors. 
#' @param X The (NT) times p design matrix, without an intercept where N=number of individuals, T=length of time series, p=number of explanatory variables.
#' @param Y The T times N panel of response where N=number of individuals, T=length of time series.
#' @param NGfactors A pre-specified number of common factors across groups (see example).
#' @param NLfactors A pre-specified number of factors in each groups (see example).
#' @param Maxit A maximum number of iterations in optimization. Default is 100.
#' @param tol Tolerance level of convergence. Default is 0.001.
#' @return A list with the following components:
#' \itemize{
#'   \item Label: The estimated group membership for each of the individuals.
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
#' @references Ando, T. and Bai, J. (2016) Panel data models with grouped factor structure under unknown group membership Journal of Applied Econometrics, 31, 163-191. 
#' @references Ando, T. and Bai, J. (2017) Clustering huge number of financial time series: A panel data approach with high-dimensional predictors and factor structures. Journal of the American Statistical Association, 112, 1182-1198. 
#' @importFrom graphics hist
#' @importFrom stats kmeans lm qnorm pnorm
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
#' PDMIFCLUST(AX,AY,R,Rs)
PDMIFCLUST <- function(X, Y, NGfactors, NLfactors, Maxit=100, tol=0.001){
  
  #Initialization
  Ngroups <- length(NLfactors)
  AY <- Y
  AX <- X
  N <- nrow(AY)
  P <- ncol(AY)
  p <- ncol(AX)
  
  Z <- AY
  VEC <- eigen(Z %*% t(Z))$vectors
  Fac <- sqrt(N)*(VEC)[,1:(sum(NLfactors))]
  L <- t(t(Fac)%*%Z/N)
  
  Km <- kmeans(L,Ngroups)
  LAB <- Km$cluster
  PP <- hist(LAB,br=0:Ngroups)$counts
  
  PredXB <- matrix(0,nrow=N,ncol=P)
  B <- matrix(0,nrow=p+1,ncol=P)
  for(j in 1:P){
    y <- Y[,j]
    X <- AX[(N*(j-1)+1):(N*j),]
    fit <- lm(y~X)
    B[,j] <- (fit$coefficients)
    PredXB[,j] <- cbind(1,X)%*%B[,j] 
  }
  
  #Local Factor
  
  Y <- AY-PredXB
  
  FS <- matrix(0,nrow=N,ncol=sum(NLfactors))
  LS <- matrix(0,nrow=P,ncol=max(NLfactors))
  PredL <- matrix(0,nrow=N,ncol=P)
  
  for(i in 1:length(NLfactors)){
    index <- subset(1:P,LAB==i)
    Z <- Y[,index]
    VEC <- eigen(Z%*%t(Z))$vectors
    Fac <- sqrt(N)*(VEC)[,1:NLfactors[i]]
    L <- t(t(Fac)%*%Z/N)
    LS[index,1:NLfactors[i]] <- L
    if(i==1){FS[,1:NLfactors[1]] <- Fac}
    if(i!=1){FS[,(sum(NLfactors[1:(i-1)])+1):(sum(NLfactors[1:i]))] <- Fac}
    PredL[,index] <- Fac%*%t(L)
  }
  
  #Global Factor
  
  Z <- AY-PredL-PredXB
  VEC <- eigen(Z%*%t(Z))$vectors; 
  Fac <- (sqrt(N)*(VEC))[,1:NGfactors]
  L <- t(t(Fac)%*%Z/N)
  PredG <- Fac%*%t(L)
  
  
  #Estimation
  
  for(ITE in 1:Maxit){
    
    LAB.old <- LAB
    B.old <- B
    
    #Beta
    
    Y <- AY-PredL-PredG
    
    for(j in 1:P){
      y <- Y[,j]
      X <- AX[(N*(j-1)+1):(N*j),]
      fit <- lm(y~X)
      B[,j] <- (fit$coefficients)
      PredXB[,j] <- cbind(1,X)%*%B[,j]
    }
    
    rm(fit)
    gc()
    
    #Lab
    
    Y <- AY-PredXB-PredG
    
    for(j in 1:P){
      
      Er <- rep(10^7,len=length(NLfactors))
      
      for(i in 1:Ngroups){
        if(NLfactors[i]!=0){
          if(i==1){Fac <- FS[,1:NLfactors[1]]}
          if(i!=1){Fac <- FS[,(sum(NLfactors[1:(i-1)])+1):(sum(NLfactors[1:i]))]}
          L <- solve(t(Fac)%*%Fac)%*%t(Fac)%*%Y[,j]
          Er[i] <- sum( (Y[,j]-Fac%*%L)^2 )
        }
        if(NLfactors[i]==0){
          Er[i] <- sum( (Y[,j])^2 )
        }
        LAB[j] <- subset(1:length(Er),Er==min(Er))
      }
      
    }
    
    #Local Factors
    
    for(i in 1:length(NLfactors)){
      if(NLfactors[i]!=0){
        index <- subset(1:P,LAB==i)
        Z <- Y[,index]
        VEC <- eigen(Z%*%t(Z))$vectors
        Fac <- sqrt(N)*(VEC)[,1:NLfactors[i]]
        L <- t(t(Fac)%*%Z/N)
        LS[index,1:NLfactors[i]] <- L
        if(i==1){FS[,1:NLfactors[1]] <- Fac}
        if(i!=1){FS[,(sum(NLfactors[1:(i-1)])+1):(sum(NLfactors[1:i]))] <- Fac}
        PredL[,index] <- Fac%*%t(L)
      }
    }
    
    #Global Factors
    
    Z <- AY-PredL-PredXB
    
    VEC <- eigen(Z%*%t(Z))$vectors; 
    Fac <- (sqrt(N)*(VEC))[,1:NGfactors]
    L <- t(Fac)%*%Z/N
    PredG <- Fac%*%L
    
    if (mean(abs(B.old - B)) <= tol) { break}
    
  }
  
  Er <- AY-PredL-PredXB-PredG
  
  Tstat  <- 0*B
  pVal  <- 0*B
  Lower05  <- 0*B
  Upper95  <- 0*B
  V <- 0*B
  B0 <- 0*B #Hypothetical
  S <- rep(0,P)
  
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
  
  
  return(list("Label"=LAB,"Coefficients"=B,"Lower05"=Lower05,"Upper95"=Upper95,"GlobalFactors"=Fac,
              "GlobalLoadings"=L,"GroupFactors"=FS,"GroupLoadings"=LS,"pval"=pVal,"Se"=V))
}


