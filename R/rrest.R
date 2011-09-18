########################################################################
#
# This file is part of my PhD Thesis
# description: reimplements the RegMCD algorithm by croux and haesbroeck
#              with L1 and L2 penalty				
# author: Moritz Gschwandtner, moritz.gschwandtner@chello.at
#
########################################################################
#require(glasso)
#require(mvoutlier)

rrest <- function(data, lambda=0.5, alpha=0.75, thresh=0.0001, maxit=10)
{
	## check requirements
	if(!is.data.frame(data) && !is.matrix(data))
		stop("parameter data has to be a matrix or data frame!")

	data <- as.data.frame(data)
		
	n <- dim(data)[1]
	h <- floor(n*alpha)
	p <- dim(data)[2]
	n_iter <- 0
	na_count <- c()
	
	if(alpha < 1) ## if h < n we have to do the whole procedure
	{
		na_count <- 0
		v_loglik <- NA
		while(is.na(v_loglik))
		{
			na_count <- na_count+1
			# use pcout for outlier detection, and focus on the inner h/2 points
			pcoutw <- pcout(data)$wfinal
			rsam <- order(pcoutw,decreasing=TRUE)[1:trunc(h/2)]
			rsub <- data[rsam,]
			## compute inverse covariance matrix using lasso penalty and mean 
			v_cov <- var(rsub)
			v_mean <- mean(rsub)
			v_gl <- glasso(v_cov, rho=lambda)
			## WARNING: glasso() sometimes produces concentration matrices with negative determinants
			## this leads to NaN in the following computation!
			## compute value of objective function
			v_loglik <- log(det(v_gl$wi)) - mean(mahalanobis(rsub, center=v_mean, cov=v_gl$wi,inverted=TRUE)) - lambda*sum(abs(v_gl$wi))
		}

		while(n_iter < maxit)
		{
			n_iter <- n_iter + 1
			## compute mahalanobis distances of all observations
			v_mahal <- mahalanobis(x=data, center=v_mean, cov=v_gl$wi, inverted=TRUE)
			v_sorted <- sort.int(v_mahal, index.return=TRUE)
			## now get the next subsample with smallest h distances
			rsam <- v_sorted$ix[1:h]
			rsub <- data[rsam,]

			## compute inverse covariance matrix using lasso penalty and mean 
			v_cov <- var(rsub)
			v_mean <- mean(rsub)
			v_gl <- glasso(v_cov, rho=lambda)
			
			## compute value of objective function and mean vector
			mds <- mahalanobis(rsub, center=v_mean, cov=v_gl$wi,inverted=TRUE)
			v_loglik_new <- log(det(v_gl$wi)) - mean(mds) - lambda*sum(abs(v_gl$wi))
		
			## this break statement assures that the algorithm does not terminate with error msg
			if(is.na( abs((v_loglik_new - v_loglik)/v_loglik) ))
				break;
				
			if(as.logical(abs((v_loglik_new - v_loglik)/v_loglik) < thresh))
				break;
				
			v_loglik <- v_loglik_new
		}
	}
	else ## alpha = 1, so glasso is just computed once
	{	
		rsam <- 1:nrow(data)
		rsub <- data[rsam,]
		v_cov <- var(rsub)
		v_mean <- mean(rsub)
		v_gl <- glasso(v_cov, rho=lambda)
		v_loglik_new <- log(det(v_gl$wi)) - mean(mahalanobis(rsub, center=v_mean, cov=v_gl$wi,inverted=TRUE)) - lambda*sum(abs(v_gl$wi))
		n_iter <- 1
		nacount <- 0
	}
	## compute consistent version
	cov_cons <- v_gl$w*median(mahalanobis(data, center=v_mean, cov=v_gl$w))/qchisq(0.5, p)
	covi_cons <- solve(cov_cons) 

	l <- list(mean=v_mean, cov=cov_cons, covi=covi_cons, cov_nocons=v_gl$w, covi_nocons=v_gl$wi, subset=rsam, loglik=v_loglik_new, niter=n_iter, nacount=na_count)
	l
}



