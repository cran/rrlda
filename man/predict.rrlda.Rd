\name{predict.rrlda}
\alias{predict.rrlda}
\alias{predict}
\title{Class Prediction for rrlda objects}
\description{
	Computes class predictions for new data based on a given rrlda object.
}
\usage{
\method{predict}{rrlda}(object, x, ...)
}
\arguments{
	\item{object}{An object of class "rrlda".}
	\item{x}{New data for which the classes are to predict}
	\item{...}{Argument used by generic function predict(object, x, ...).}
}
\details{
	Based on the estimated inverse covariance matrix and the mean of each group (stored in object), discriminant values are computed. An observations is classified as group k, if the corresponding discriminant value is a minimum.  
}
\value{
	\item{class}{Class prediction for each observation.}
	\item{posterior}{Discriminant values.}
}
\examples{
	data(iris)
	x <- iris[,1:4]
	rr <- rrlda(x, grouping=as.numeric(iris[,5]), lambda=0.2, alpha=0.75) ## perform rrlda
	pred <- predict(rr, x) ## predict 
	table(as.numeric(pred$class), as.numeric(iris[,5])) ## show errors
}
\keyword{file}

