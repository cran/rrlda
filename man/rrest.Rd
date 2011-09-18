\name{rrest}
\alias{rrest}
\title{Robust Regularized Estimator (RegMCD) for location and inverse scatter}
\description{
  Computes the Robust Regularized Estimator for location and inverse scatter. 
}
\usage{
rrest(data, lambda, alpha, thresh, maxit)
}
\arguments{
	\item{data}{Matrix or data.frame of observations}
	\item{lambda}{Penalty parameter which controls the sparseness of the resulting inverse scatter matrix. Default is 0.5}
	\item{alpha}{Robustness parameter which specifies the amount of observations to be included in the computations. Default is 0.75}
	\item{thresh}{Threshold value controlling the convergence of the iterative algorithm. Default is 0.0001. In most cases this argument does not have to be supplied.}
	\item{maxit}{Maximum number of iterations of the algorithm. Default is 10.}
}
\details{
	The Robust Regularized Estimator computes a sparse inverse covariance matrix of the given observations by maximization of a penalized likelihood function. The sparseness is controlled by a penalty parameter lambda. Possible outliers are dealt with by a robustness parameter alpha which specifies the amount of observations for which the likelihood function is maximized. 
}
\value{

	\item{mean}{The resulting location estimate.}
	\item{cov}{The resulting covariance estimate.}
	\item{covi}{The resulting inverse covariance estimate.}
	\item{subset}{An index vector specifying the data subset used (see robustness parameter alpha).}
	\item{loglik}{The maximized (log-)likelihood value.}
	\item{niter}{The number of iterations}
}
\examples{
	x <- cbind(rnorm(100), rnorm(100), rnorm(100)) # use first group only
	rr <- rrest(x, lambda=0.2, alpha=0.75)
	rr$cov ## estimated cov matrix
}
\keyword{file}

