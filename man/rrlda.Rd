\name{rrlda}
\alias{rrlda}
\title{Robust Regularized Linear Discriminant Analysis}
\description{
	Performs Robust Regularized Linear Discriminant Analysis.
}
\usage{
rrlda(x, grouping, prior=NULL, lambda=0.5, hp=0.75, nssamples=30, maxit=50, penalty="L2")
}
\arguments{
	\item{x}{Matrix or data.frame of observations.}
	\item{grouping}{Grouping variable. A vector of numeric values >= 1 is recommended. Length has to correspond to nrow(x).}
	\item{prior}{Vector of prior probabilities for each group. If not supplied the priors are computed from the data.}
	\item{lambda}{Penalty parameter which controls the sparseness of the resulting inverse scatter matrix. Default is 0.5}
	\item{hp}{Robustness parameter which specifies the amount of observations to be included in the computations. Default is 0.75}
	\item{nssamples}{Number of start samples to be user for iterated estimations.}
	\item{maxit}{Maximum number of iterations of the algorithm. Default is 10.}
	\item{penalty}{Type of penalty to be applied. Possible values are "L1" and "L2".}
}
\details{
	Performs Robust Regularized Discriminant Analysis using a sparse estimation of the inverse covariance matrix. The sparseness is controlled by a penalty parameter lambda. Possible outliers are dealt with by a robustness parameter alpha which specifies the amount of observations for which the likelihood function is maximized. 
}
\value{
	An object of class "rrlda" is returned which can be used for class prediction (see predict()).
prior=prior, counts=counts, means=means, cov=covm, covi=covi, lev=lev, n=n, h=h, bic=bic, loglik=loglik, nonnuls=nonnuls, subs=est$subset

	\item{prior}{Vector of prior probabilities.}
	\item{counts}{Number of obervations for each group.}
	\item{means}{Estimated mean vectors for each group.}
	\item{covi}{Estimated (common) inverse covariance matrix.}
	\item{lev}{Levels. Corresponds to the groups.}
	\item{n}{Number of observations.}
	\item{h}{Number of observations included in the computations (see robustness parameter alpha).}
	\item{bic}{Adapted bic value. Can be used for optimal selection of lambda}
	\item{loglik}{The maximized (log-)likelihood value.}
	\item{df}{Degrees of freedom of the estimated inverse covariance matrix.}
	\item{subs}{An index vector specifying the data subset used (see robustness parameter alpha).}
}
\examples{
	data(iris)
	x <- iris[,1:4]
	rr <- rrlda(x, grouping=as.numeric(iris[,5]), lambda=0.2, hp=0.75) ## perform rrlda
	pred <- predict(rr, x) ## predict 
	table(as.numeric(pred$class), as.numeric(iris[,5])) ## show errors
}
\keyword{file}

