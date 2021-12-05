#' HYPTEST
#'
#' This function undergoes hypothesis testing for regression coefficients obtained from the various functions in the package.
#' @param B A dataframe of Coefficients as obtained in the output of any function in the package.
#' @param B0 A dataframe of hypothetical coefficients to be evaluated in the test. (nrows should match number of variables and ncols should match number of individuals) 
#' @param Se A dataframe of Standard Errors as obtained in the output of any function in the package.
#' @param test A string to determine what kind of test to run ("two" for two-tailed, "right" for right-tailed and "left for left-tailed).
#' @param individuals A list of individuals whose coefficients are to be tested. Default is all individuals in the B dataframe.
#' @param variables A list of variables whose coefficients are to be tested. Default is all variables in the B dataframe.
#' @return A dataframe of p-values resulting from each individual test.
#' @importFrom stats pnorm
#' @export
#' @examples 
#' N <- 400
#' P <- 300
#' 
#' R <- 3
#' p <- 3
#' 
#' LAM <- matrix(rnorm(P*R,0,1),nrow=P,ncol=R)
#' FAC <- matrix(rnorm(N*R,0,1),nrow=N,ncol=R)
#' AY <- FAC%*%t(LAM)
#' 
#' AX <- matrix(rnorm(p*P*N,0,1),nrow=P*N)
#' for(j in 1:P){
#'   AX[(N*(j-1)+1):(N*j),1] <- AX[(N*(j-1)+1):(N*j),1]+0.15*AY[,j]
#'   AX[(N*(j-1)+1):(N*j),3] <- AX[(N*(j-1)+1):(N*j),3]-0.21*AY[,j]
#' }
#' 
#' B <- c(-1,2,-1)
#' AB <- B%*%t(rep(1,len=P))+matrix(0.1*runif(p*P,-1,1),p,P)
#' for(j in 1:P){
#'   AY[,j] <- AY[,j]+(AX[(N*(j-1)+1):(N*j),])%*%AB[,j]
#' }
#' 
#' PROB <- exp(AY)/(1+exp(AY))
#' AY <- trunc(matrix(runif(N*P,0,1),ncol=P)+PROB)
#' fit <- PDMIFLOGIT(AX,AY,R)
#' HYPTEST(fit$Coefficients,data.frame(c(0,1),c(-1,2)),fit$Se,"two",c(1,2),c(1,3))
HYPTEST <- function(B,B0,Se,test="two",individuals=seq(1,ncol(B)),variables=seq(1,nrow(B))){

  if (ncol(B0) != length(individuals)){
  stop("ncols of hypothetical beta dataframe and length of individuals vector do not match")
  }
  if (nrow(B0) != length(variables)){
    stop("nrows of hypothetical beta dataframe and length of variables vector do not match")
  }
  B<-data.frame(B[variables,individuals])
  Se<-data.frame(Se[variables,individuals])
  Tstat <- (B-B0)/Se
  if(test=="two"){
  pVal <- 2*apply(-abs(Tstat),2,pnorm)
  }
  else if(test=="right"){
    pVal <- apply(-Tstat,2,pnorm)
  }
  else{pVal <- apply(Tstat,2,pnorm)}
  return(pVal)
}
