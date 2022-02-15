#' PDMIFCLUSTGLM
#'
#' Under a pre-specified number of groups and the number of common factors, this function implements clustering for N individual units by nonlinear heterogeneous panel data models with interactive effects. 
#' Exponential family of distributions are used
#' Each of individuals in the group are subject to the group-specific unobserved common factors. 
#' @param X The (NT) times p design matrix, without an intercept where N=number of individuals, T=length of time series, p=number of explanatory variables.
#' @param Y The T times N panel of response where N=number of individuals, T=length of time series.
#' @param FAMILY A description of the error distribution and link function to be used in the model just like in glm functions.
#' @param NLfactors A pre-specified number of factors in each groups (see example).
#' @param Maxit A maximum number of iterations in optimization. Default is 100.
#' @param tol Tolerance level of convergence. Default is 0.001.
#' @return A list with the following components:
#' \itemize{
#'   \item Label: The estimated group membership for each of the individuals.
#'   \item Coefficients: The estimated heterogeneous coefficients.
#'   \item Lower05: Lower end (5%) of the 90% confidence interval of the regression coefficients.
#'   \item Upper95: Upper end (95%) of the 90% confidence interval of the regression coefficients.
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
#' fit <- PDMIFCLUSTGLM(data6X,data6Y,binomial(link=logit),c(1,1,1),3,0.5)
PDMIFCLUSTGLM <- function(X, Y, FAMILY, NLfactors, Maxit=100, tol=0.001){
  set.seed(1221)
  N <- nrow(Y) #length of time series
  p <- ncol(X) #number of explanatory variables
  P <- ncol(Y) #total number of individuals
  AY<-Y
  AX<-X

  Km <- kmeans(t(AY),length(NLfactors))
  LAB <- Km$cluster

  #Initialization
  
  PredXB <- matrix(0,nrow=N,ncol=P)
  B <- matrix(0,nrow=p,ncol=P)
  
  for(j in 1:P){
    y <- AY[,j]
    X <- AX[(N*(j-1)+1):(N*j),]
    fit <- glm(y~X+0,family=binomial(link=logit))
    B[,j] <- (fit$coefficients)
    PredXB[,j] <- X%*%B[,j]
  }
  
  FS <- matrix(0,nrow=N,ncol=sum(NLfactors))
  LS <- matrix(0,nrow=P,ncol=max(NLfactors))
  PredL <- matrix(0,nrow=N,ncol=P)
  
  for(i in 1:length(NLfactors)){
    index <- subset(1:(P),LAB==i)
    Z <- (AY-PredXB)[,index]
    NG <- length(index)
    VEC <- eigen(Z%*%t(Z))$vectors
    Ftemp <- sqrt(NG)*(VEC)[,1:NLfactors[i]]
    Ltemp <- t(t(Ftemp)%*%Z/NG)
    LS[index,1:NLfactors[i]] <- Ltemp
    if(i==1){FS[,1:NLfactors[1]] <- Ftemp}
    if(i!=1){FS[,(sum(NLfactors[1:(i-1)])+1):(sum(NLfactors[1:i]))] <- Ftemp}
    PredL[,index] <- Ftemp%*%t(Ltemp)
  }
  
  AB <- rbind(B,t(LS))
  B.old <- B
  
  for(j in 1:P){
    i <- LAB[j]
    if(i==1){Ftemp <- FS[,1:NLfactors[1]]}
    if(i!=1){Ftemp <- FS[,(sum(NLfactors[1:(i-1)])+1):(sum(NLfactors[1:i]))]}
    y <- AY[,j]
    X <- cbind(AX[(N*(j-1)+1):(N*j),],Ftemp)
    fit <- glm(y~X+0,family=binomial(link=logit))
    AB[1:(p+NLfactors[i]),j] <- fit$coefficients
    LS[j,1:NLfactors[i]] <- AB[(p+1):(p+NLfactors[i]),j]
    B[,j] <- AB[1:p,j]
    PredXB[,j] <- AX[(N*(j-1)+1):(N*j),]%*%AB[1:p,j]
  }
  
  for(ITE in 1:Maxit){
    
    for(i in 1:length(NLfactors)){
      
      index <- subset(1:(P),LAB==i)
      Ltemp <- LS[index,1:NLfactors[i]]
      
      Ftemp <- matrix(0,N,NLfactors)
      for(k in 1:N){
        y <- AY[k,index]
        X <- Ltemp
        fit <- glm(y~X+offset(PredXB[k,index])+0,family=binomial(link=logit))
        Ftemp[k,] <- fit$coefficients
      }
      
      NG <- sum(LAB==length(NLfactors))
      QRF <- qr(Ftemp)
      Ftemp <- sqrt(NG)*qr.Q(QRF)
      if(i==1){FS[,1:NLfactors[1]] <- Ftemp}
      if(i!=1){FS[,(sum(NLfactors[1:(i-1)])+1):(sum(NLfactors[1:i]))] <- Ftemp}
      
    }
    
    for(j in 1:P){
      i <- LAB[j]
      if(i==1){Ftemp <- FS[,1:NLfactors[1]]}
      if(i!=1){Ftemp <- FS[,(sum(NLfactors[1:(i-1)])+1):(sum(NLfactors[1:i]))]}
      y <- AY[,j]
      X <- cbind(AX[(N*(j-1)+1):(N*j),],Ftemp)
      fit <- glm(y~X+0,family=binomial(link=logit))
      AB[1:(p+NLfactors[i]),j] <- fit$coefficients
      LS[j,1:NLfactors[i]] <- AB[(p+1):(p+NLfactors[i]),j]
      B[,j] <- AB[1:p,j]
      PredXB[,j] <- AX[(N*(j-1)+1):(N*j),]%*%AB[1:p,j]
    }
  
    for(j in 1:P){
      Er <- rep(10^7,len=length(NLfactors))
      for(i in 1:length(NLfactors)){
        if(i==1){Ftemp <- FS[,1:NLfactors[1]]}
        if(i!=1){Ftemp <- FS[,(sum(NLfactors[1:(i-1)])+1):(sum(NLfactors[1:i]))]}
        y <- AY[,j]; X <- Ftemp
        fit <- glm(y~X+offset(PredXB[,j])+0,family=binomial(link=logit))
        Er[i] <- fit$aic
      }
      LAB[j] <- subset(1:length(Er),Er==min(Er))
    }
    if(mean(abs(B.old-B))<=tol){break}
    B.old <- B

  }#ITE
  


Er <- AY-PredL-PredXB

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
  V[,i] <- sqrt(diag(S[i]*t(X)%*%X))
  Lower05[,i]  <- B[,i]+qnorm(0.05)*V[,i]/sqrt(N)
  Upper95[,i]  <- B[,i]+qnorm(0.95)*V[,i]/sqrt(N)
  Tstat[,i] <- sqrt(N)*(B[,i]-B0[,i])/V[,i]
  pVal[,i] <- 2*pnorm(-abs(Tstat[,i]))
}
cat("Call:
PDMIFCLUSTGLM(X, Y, FAMILY =",FAMILY$family,FAMILY$link,", NLfactors =",NLfactors,", Maxit =",Maxit,", tol =",tol,")
  
N =",P,", T =",N,", p =",p,"

Fit includes coefficients, confidence interval, group factors, group loadings, 
p-values and standard errors.
")

invisible(list("Label"=LAB,"Coefficients"=B,"Lower05"=Lower05,"Upper95"=Upper95,
               "GroupFactors"=FS,"GroupLoadings"=LS,"pval"=pVal,"Se"=V))
}
