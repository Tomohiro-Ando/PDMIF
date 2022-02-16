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
#' fit <- PDMIFCLUST(data5X,data5Y,2,c(2,2,2),20,0.5)
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
  Ftemp <- sqrt(N)*(VEC)[,1:(sum(NLfactors))]
  Ltemp <- t(t(Ftemp)%*%Z/N)
  
  Km <- kmeans(Ltemp,Ngroups)
  LAB <- Km$cluster
  PP <- hist(LAB,br=0:Ngroups,plot=FALSE)$counts
  
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
  
  for(i in 1:Ngroups){
    index <- subset(1:P,LAB==i)
    Z <- Y[,index]
    VEC <- eigen(Z%*%t(Z))$vectors
    Ftemp <- sqrt(N)*(VEC)[,1:NLfactors[i]]
    Ltemp <- t(t(Ftemp)%*%Z/N)
    LS[index,1:NLfactors[i]] <- Ltemp
    if(i==1){FS[,1:NLfactors[1]] <- Ftemp}
    if(i!=1){FS[,(sum(NLfactors[1:(i-1)])+1):(sum(NLfactors[1:i]))] <- Ftemp}
    PredL[,index] <- Ftemp%*%t(Ltemp)
  }
  
  #Global Factor
  
  Z <- AY-PredL-PredXB
  VEC <- eigen(Z%*%t(Z))$vectors; 
  Ftemp <- (sqrt(N)*(VEC))[,1:NGfactors]
  Ltemp <- t(t(Ftemp)%*%Z/N)
  FG <- Ftemp
  LG <- Ltemp
  PredG <- FG%*%t(LG)
  
  
  #Estimation
  
  for(ITE in 1:Maxit){

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
      
      Er <- rep(10^7,len=Ngroups)
      
      for(i in 1:Ngroups){
        if(NLfactors[i]!=0){
          if(i==1){Ftemp <- FS[,1:NLfactors[1]]}
          if(i!=1){Ftemp <- FS[,(sum(NLfactors[1:(i-1)])+1):(sum(NLfactors[1:i]))]}
          Ltemp <- solve(t(Ftemp)%*%Ftemp)%*%t(Ftemp)%*%Y[,j]
          Er[i] <- sum( (Y[,j]-Ftemp%*%Ltemp)^2 )
        }
        if(NLfactors[i]==0){
          Er[i] <- sum( (Y[,j])^2 )
        }
        LAB[j] <- subset(1:length(Er),Er==min(Er))
      }
      
    }
    
    #Local Factors
    
    for(i in 1:Ngroups){
      if(NLfactors[i]!=0){
        index <- subset(1:P,LAB==i)
        Z <- Y[,index]
        VEC <- eigen(Z%*%t(Z))$vectors
        Ftemp <- sqrt(N)*(VEC)[,1:NLfactors[i]]
        Ltemp <- t(t(Ftemp)%*%Z/N)
        LS[index,1:NLfactors[i]] <- Ltemp
        if(i==1){FS[,1:NLfactors[1]] <- Ftemp}
        if(i!=1){FS[,(sum(NLfactors[1:(i-1)])+1):(sum(NLfactors[1:i]))] <- Ftemp}
        PredL[,index] <- Ftemp%*%t(Ltemp)
      }
    }
    
    #Global Factors
    
    Z <- AY-PredL-PredXB
    
    VEC <- eigen(Z%*%t(Z))$vectors; 
    Ftemp <- (sqrt(N)*(VEC))[,1:NGfactors]
    Ltemp <- t(t(Ftemp) %*% Z/N)
    FG <- Ftemp
    LG <- Ltemp
    PredG <- FG%*%t(LG)
    
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
    Lower05[,i]  <- B[,i]+qnorm(0.05)*V[,i]/sqrt(N)
    Upper95[,i]  <- B[,i]+qnorm(0.95)*V[,i]/sqrt(N)
    Tstat[,i] <- sqrt(N)*(B[,i]-B0[,i])/V[,i]
    pVal[,i] <- 2*pnorm(-abs(Tstat[,i]))
  }
  message("Call:
PDMIFCLUST(X, Y, NGfactors =",NGfactors,", NLfactors =",NLfactors,", Maxit =",Maxit,", tol =",tol,")
  
N =",P,", T =",N,", p =",p,"

Fit includes coefficients, confidence interval, global factors, global loadings,
    group factors, group loadings, p-values and standard errors.
")
  
  invisible(list("Label"=LAB,"Coefficients"=B,"Lower05"=Lower05,"Upper95"=Upper95,"GlobalFactors"=FG,
              "GlobalLoadings"=LG,"GroupFactors"=FS,"GroupLoadings"=LS,"pval"=pVal,"Se"=V))
}


