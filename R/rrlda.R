########################################################################
#
# This file is part of my PhD Thesis
# description: It contains the functions rrlda and predict.rrlda 
#              performing robust regularized linear discriminant 
#              analysis and the prediction of new observations based on
#              a model of class rrlda
# author: Moritz Gschwandtner, moritz.gschwandtner@chello.at
#
########################################################################

#source("rrest.R")
#require(pcaPP)
#source("haesbroeck.R")

rrlda <- function(x, grouping, prior=NULL, lambda=0.5, alpha=0.75, maxit=50)
{
	## check requirements
    if (is.null(dim(x))) 
	{
        stop("'x' is not a matrix or data.frame")
	}
	x <- as.matrix(x)
    if (any(!is.finite(x))) 
	{
        stop("Infinite, NA or NaN values in 'x'")
	}
	## get dimensions
	n <- nrow(x)
    p <- ncol(x)
	h <- floor(n*alpha)
	## check if length(x) == length(grouping)
    if (n != length(grouping)) 
	{
        stop("nrow(x) and length(grouping) are different")	
	}
	## is grouping a factor?
	if(!is.factor(grouping))
	{
		grouping <- as.factor(grouping)
	}
	## get levels, counts ands proportions
	lev = levels(grouping)
	k = nlevels(grouping)
	counts = table(grouping)
	proportions <- counts/n
	
	## check prior probabilities
	if (!missing(prior)) 
	{
        if (any(prior < 0) || round(sum(prior), 5) != 1) 
		{
            stop("invalid prior")
		}
        if (length(prior) != k) 
		{
            stop("'prior' is of incorrect length")
		}
    }
	else
	{
		prior <- proportions
	}
	## compute common covariance matrix
	covm <- matrix(rep(0, p*p), ncol=p)
	## non consistent version used for bic computation
	covm_nocons <- matrix(rep(0, p*p), ncol=p)
	## subset of all observations used  
	subs <- c()
	## corresponding subset of grouping variable
	gr_subs <- c()

	## compute location estimates based on RegMCD
	m <- t(x)
	for(i in 1:k)
	{
		# we could either calculate the RegMCD estimator for each group (takes very long!) or simply center data by the median
		#est <- rrest(x[grouping==lev[i],], lambda=lambda, alpha=alpha, penalty.fct=penalty, nssamples=nssamples, maxit=maxit)
		#m[,grouping==lev[i]] <- est$mean
		#m[,grouping==lev[i]] <- apply(x[grouping==lev[i],], 2, median)
		## L1 median: (WARNING: needs at least 2 observations in each group! thats the reason for the if)
		if(!is.vector(x[grouping==lev[i],]))
		{
			m[,grouping==lev[i]] <- l1median_NLM(x[grouping==lev[i],])$par
		}
		else
		{
			m[,grouping==lev[i]] <- as.vector(x[grouping==lev[i],])
		}
	}
	m <- t(m)
	z <- x - m
	
	## compute covariance matrix of centered observations
	est <- rrest(z, lambda=lambda, alpha=alpha, maxit=maxit)
	subs <- x[est$subset,]
	gr_subs <- grouping[est$subset]
	covm <- est$cov
	covm_nocons <- est$cov_nocons
	
	# compute final group means
	m_final <- t(t(m) + est$mean)
	
	# create index vector with indices of first, second, etc. level; purpose: sort means afterwards
	ind_vector <- c()
	for(i in 1:length(lev))
	{
		ind_vector <- c(ind_vector, which(grouping == lev[i]))
	}
	m_final <- m_final[ind_vector, ]
	
	# now remove double entries, so we have just one row per group
	# rows are ordered according to levels; in general this isn't the same order as the groups occur in x
	means <- unique(m_final)
	## compute concentration matrix
	rownames(means) <- lev
	covi <- solve(covm)
	covi_nocons <- solve(covm_nocons)
	
	## compute bic without penalty term in likelihood function
	bic <- h*log(det(covi_nocons))
	for(i in 1:k)
	{
		bic <- bic - sum(mahalanobis(subs[gr_subs==lev[i],], center=means[i,], cov=covi_nocons, inverted=TRUE))
	}
	nonnuls <- sum(covi_nocons[upper.tri(covi_nocons,diag=TRUE)]!=0)
	loglik <- (-1)*bic
	bic <- loglik + log(h)*(k*p + nonnuls)
	
	cl <- match.call()
    cl[[1L]] <- as.name("rrlda")
	
	ret = structure(list(call=cl, prior=prior, counts=counts, means=means, cov=covm, covi=covi, lev=lev, n=n, h=h, bic=bic, loglik=loglik, nonnuls=nonnuls, subs=est$subset), class="rrlda")
	ret
}

predict.rrlda <- function(object, x, ...)
{
	## check requirements
	if(class(object) != "rrlda")
		stop("object has to be of type rrlda")
		
	## get dimension	
	p = ncol(object$means)

	if(!is.vector(x) && !is.matrix(x) && !is.data.frame(x))
		stop("'x' has to be a vector, matrix or data.frame")
	
	n <- 0
	if(is.vector(x))
	{
		if(length(x) != p)
			stop("'x' has to be of length ", p)

		y <- NULL
		for(i in 1:length(x))
			y <- cbind(y, x[i])
		x <- y
		n <- 1
	}
	else
	{
		if(ncol(x) != p)
			stop("'x' has to be of dimension ", p)
			
		x <- as.matrix(x)
		n <- nrow(x)
	}
	
	#prior_term <- -2*log(t(as.matrix(object$prior)))
	#mean_term <- diag(object$means %*% object$covi %*% t(object$means))
	
	## create matrices for calculation
	#for(i in 1:n)
	#{
	#	if(i != 1)
	#	{
	#		prior_term <- rbind(prior_term, -2*log(t(as.matrix(object$prior))))
	#		mean_term <- rbind(mean_term, diag(object$means %*% object$covi %*% t(object$means)))
	#	}
	#}
	
	#probs <- -2*t(object$means %*% object$covi %*% t(x)) + mean_term + prior_term
	
	probs <- t(-2*object$means %*% object$covi %*% t(x) + diag(object$means %*% object$covi %*% t(object$means)) - 2*log(object$prior[[1]]))
	cl = object$lev[apply(probs, 1, which.min)]
	
	ret = list(class=cl, posterior=probs)
	ret
}
