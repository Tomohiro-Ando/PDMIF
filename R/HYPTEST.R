#' HYPTEST
#'
#' This function undergoes hypothesis testing for regression coefficients obtained from the various functions in the package.
#' @param B A dataframe of Coefficients as obtained in the output of any function in the package.
#' @param B0 A dataframe of hypothetical coefficients to be evaluated in the test. (nrows should match number of variables and ncols should match number of individuals) 
#' @param Se A dataframe of Standard Errors as obtained in the output of any function in the package.
#' @param test A string to determine what kind of test to run ("two" for two-tailed, "right" for right-tailed and "left for left-tailed).
#' @param variables A list of variables whose coefficients are to be tested. Default is all variables in the B dataframe.
#' @param individuals A list of individuals whose coefficients are to be tested. Default is all individuals in the B dataframe.
#' @return A dataframe of p-values resulting from each individual test.
#' @importFrom stats pnorm
#' @export
#' @examples 
#' fit <- PDMIFLOGIT(data2X,data2Y,2)
#' HYPTEST(fit$Coefficients,data.frame(c(0,1),c(-1,2)),fit$Se,"two",c(1,3),c(1,2))
HYPTEST <- function(B,B0,Se,test="two",variables=seq(1,nrow(B)),individuals=seq(1,ncol(B))){

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
  pVal=data.frame(pVal)
  colnames(pVal)=paste("i",individuals)
  rownames(pVal)=paste("v",variables)
  return(pVal)
}
