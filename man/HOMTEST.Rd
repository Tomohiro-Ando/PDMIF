% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HOMTEST.R
\name{HOMTEST}
\alias{HOMTEST}
\title{HOMTEST}
\usage{
HOMTEST(X, Y, Nfactors, Maxit = 100, tol = 0.001)
}
\arguments{
\item{X}{The (NT) times p design matrix, without an intercept where N=number of individuals, T=length of time series, p=number of explanatory variables.}

\item{Y}{The T times N panel of response where N=number of individuals, T=length of time series.}

\item{Nfactors}{A pre-specified number of common factors.}

\item{Maxit}{A maximum number of iterations in optimization. Default is 100.}

\item{tol}{Tolerance level of convergence. Default is 0.001.}
}
\value{
A list with the following components:
\itemize{
\item Coefficients: The estimated heterogeneous coefficients.
\item Factors: The estimated common factors across groups.
\item Loadings: The estimated factor loadings for the common factors.
\item pvalue: The p-value of the homogeneity test.
}
}
\description{
This function tests homogeneity of the regression coefficients in heterogeneous panel data models with interactive effects.
}
\examples{
fit <- HOMTEST(data1X,data1Y,2,20,0.5)
}
\references{
Ando, T. and Bai, J. (2015) A simple new test for slope homogeneity in panel data models with interactive effects. Economics Letters, 136, 112-117.
}
