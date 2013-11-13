#
#	Implementation of Polya-Aeppli distribution (also known as geometric-Poisson)
#
#	dPolyaAeppli(x, lambda, prob, log = FALSE)
#	pPolyaAeppli(q, lambda, prob, lower.tail = TRUE, log.p = FALSE)
#	qPolyaAeppli(p, lambda, prob, lower.tail = TRUE, log.p = FALSE)
#	rPolyaAeppli(n, lambda, prob)
#
#	Conrad Burden	November, 2013
#
########################################################################################
#
#	Probability density of the Polya-Aeppli distribution
#
	dPolyaAeppli <- function(x, lambda, prob, log=FALSE){
		if(!is.numeric(x)){
			stop("Non-numeric argument to mathematical function \n")
			} 
		if(!is.numeric(lambda)){
			stop("Non-numeric parameter lambda \n")
			} 
		if(lambda <= 0){
			stop("parameter lambda must be > 0 \n")
			} 
		if(!is.numeric(prob)){
			stop("Non-numeric parameter prob \n")
			} 
		if(prob < 0 | prob >=1){
			stop("parameter prob must be between 0 and 1 \n")
			}
		if(all(x == Inf)){
			if(log){
				return(rep(-Inf, length(x)))
				}else{
				return(rep(0, length(x)))
				}
			} 
		xMax <- ceiling(max(x[x!=Inf]))
		isValid <- (x >= 0) & (is.wholenumber(x)) & (x < Inf)
		if(log){
			dPolyaAeppli <- rep(-Inf, length(x))
			lPArray <- lPolyaAeppliArray(xMax, lambda, prob)
			dPolyaAeppli[isValid] <- lPArray[x[isValid] + 1]
			}else{
			dPolyaAeppli <- rep(0, length(x))
			dPArray <- exp(lPolyaAeppliArray(xMax, lambda, prob))
			dPolyaAeppli[isValid] <- dPArray[x[isValid] + 1]
			}
#
		warningNeeded <- !isValid & x!=Inf & x!=-Inf
		if(length(x[warningNeeded]) > 0){
			warning(paste("\n non-positive-integer x = ", x[warningNeeded], sep=""))
			}
#
		return(dPolyaAeppli)		
		}
#
########################################################################################
#
#	Cumulative probability function of the Polya-Aeppli distribution
#
	pPolyaAeppli <- function(q, lambda, prob, lower.tail=TRUE, log.p=FALSE){
		x <- q
		if(!is.numeric(x)){
			stop("Non-numeric argument to mathematical function \n")
			} 
		if(!is.numeric(lambda)){
			stop("Non-numeric parameter lambda \n")
			} 
		if(lambda <= 0){
			stop("parameter lambda must be > 0 \n")
			} 
		if(!is.numeric(prob)){
			stop("Non-numeric parameter prob \n")
			} 
		if(prob < 0 | prob >=1){
			stop("parameter prob must be between 0 and 1 \n")
			} 
		if(all(x == Inf)){
			if( log.p &  lower.tail){return(rep(0, length(x)))}
			if(!log.p &  lower.tail){return(rep(1, length(x)))}
			if( log.p & !lower.tail){return(rep(-Inf, length(x)))}
			if(!log.p & !lower.tail){return(rep(0, length(x)))}
			} 
		xFloor <- floor(x)
		xMax <- max(xFloor[x!=Inf])
		PAmean <- lambda/(1 - prob)
#
		if(lower.tail){
			lArray <- lPolyaAeppliArray(xMax, lambda, prob)
			gg <- gArray(lArray)
			lPA <- rep(-Inf,length(x))
			lPA[xFloor>=0 & x!=Inf] <- gg[xFloor[xFloor>=0 & x!=Inf] + 1]
			lPA[x==Inf] <- 0
			}
		if(!lower.tail & xMax>PAmean){
			lArray <- lPolyaAeppliArray(xMax, lambda, prob)
			hTop <- logTailPA(xMax, lambda, prob)
			hh <- hArray(hTop, lArray)
			lPA <- rep(0, length(x))
			lPA[xFloor>=0 & x!=Inf] <- hh[xFloor[xFloor>=0 & x!=Inf] + 1]
			lPA[x==Inf] <- -Inf			
			}
		if(!lower.tail & xMax<=PAmean){
			lPA <- log1p(-pPolyaAeppli(x, lambda, prob))
			}
		if(!log.p){
			pPolyaAeppli <- exp(lPA)
			}else{
			pPolyaAeppli <- lPA
			}
		return(pPolyaAeppli)		
		}
#
#######################################################################################
#
#	Quantile function of the Polya-Aeppli distribution
#
	qPolyaAeppli <- function(p, lambda, prob, lower.tail=TRUE, log.p=FALSE){
		if(!is.numeric(p)){
			stop("Non-numeric argument to mathematical function \n")
			} 
		if(!is.numeric(lambda)){
			stop("Non-numeric parameter lambda \n")
			} 
		if(lambda <= 0){
			stop("parameter lambda must be > 0 \n")
			} 
		if(!is.numeric(prob)){
			stop("Non-numeric parameter prob \n")
			} 
		if(prob < 0 | prob >=1){
			stop("parameter prob must be between 0 and 1 \n")
			} 
#
		qPolyaAeppli <- array(dim=length(p))
		needsCalculating <- rep(FALSE, length(p)) 
#
		if( log.p &  lower.tail){
			qPolyaAeppli[p==-Inf] <- 0
			qPolyaAeppli[p<=0 & p>=(- 10*.Machine$double.eps)] <- Inf
			qPolyaAeppli[p>0] <- NaN
			needsCalculating[p<(- 10*.Machine$double.eps) & p>-Inf] <- TRUE
			}
		if(!log.p &  lower.tail){
			qPolyaAeppli[p==0] <- 0
			qPolyaAeppli[p<=1 & p>=(1 - 10*.Machine$double.eps)] <- Inf
			qPolyaAeppli[p<0 | p>1] <- NaN
			needsCalculating[p>0 & p<(1 - 10*.Machine$double.eps)] <- TRUE
			}
		if( log.p & !lower.tail){
			qPolyaAeppli[p==0] <- 0
			qPolyaAeppli[p==-Inf] <- Inf
			qPolyaAeppli[p>0] <- NaN
			needsCalculating[p<0 & p>-Inf] <- TRUE
			}
		if(!log.p & !lower.tail){
			qPolyaAeppli[p==1] <- 0
			qPolyaAeppli[p>=0 & p<=(10*.Machine$double.eps)] <- Inf
			qPolyaAeppli[p<(10*.Machine$double.eps) | p>1] <- NaN
			needsCalculating[p>0 & p<1] <- TRUE
			}
#		
		if(length(p[needsCalculating])>0){
			if(log.p){
				logP <- p[needsCalculating]
				}else{
				logP <- log(p[needsCalculating])
				}
#
#		Upper bound qmax on highest (non-infinite) quantile: use quantile 
#			of gamma distribution plus 1 standard deviation
#
			if(lower.tail){
				logPExtreme <- max(logP)
				}else{
				logPExtreme <- min(logP)
				}
#
			PAmean <- lambda/(1 - prob)
			PAvar <- PAmean*(1 + prob)/(1 - prob)
			scale <- PAvar/PAmean
			shape <- PAmean/scale
			qMax <- ceiling(qgamma(logPExtreme, shape=shape, scale=scale, 
						log.p=TRUE, lower.tail=lower.tail) + sqrt(PAvar))
#
			pArray <- pPolyaAeppli(0:qMax, lambda, prob, 
										log.p=TRUE, lower.tail=lower.tail) 
#
			qPA <- c()
			if(lower.tail){
				for (thisP in logP){
					qPA <- c(qPA, max(which((pArray + .Machine$double.eps) < thisP), 0))
					}
				}else{
				for (thisP in logP){
					qPA <- c(qPA, max(which((pArray - .Machine$double.eps) > thisP), 0))
					}
				}
			qPolyaAeppli[needsCalculating] <- qPA
			}
#
		if(any(is.nan(qPolyaAeppli))){
			warning(paste("NaNs produced \n  ", sep=""))
			}
#
		return(qPolyaAeppli)
		}
#
########################################################################################
#
#	Generate a random sample from the Polya-Aeppli distribution
#
	rPolyaAeppli <- function(n, lambda, prob){
		if(!is.numeric(n)){
			stop("Non-numeric argument to mathematical function \n")
			} 
		if(!is.numeric(lambda)){
			stop("Non-numeric parameter lambda \n")
			} 
		if(lambda <= 0){
			stop("parameter lambda must be > 0 \n")
			} 
		if(!is.numeric(prob)){
			stop("Non-numeric parameter prob \n")
			} 
		if(prob < 0 | prob >=1){
			stop("parameter prob must be between 0 and 1 \n")
			} 
		rPolyaAeppli <- array(dim=n)
		for(i in 1:n){
			nSample <- rpois(1, lambda)
			rPolyaAeppli[i] <- sum(rgeom(nSample, (1 - prob))) + nSample
			}
			rPolyaAeppli
		}
#
#######################################################################################
#
#	Array whose values are the log of the probability function, log(Pr(X=x)), 
#		for the Polya-Aeppli distribution for x from 0 up to max(2, xMax).   
#		(note that the array index is 1 more than x, i.e. lArray[x] = Prob(X=(x-1))) 
#	Uses an iterative formula based on Eq.(9.165) from p379, Johnson, Kotz and Kemp, 
#		also used in Nuel, J. Stat. Comp. Sim. 78 (2008) 385-394
#
	lPolyaAeppliArray <- function(xMax, lambda, prob){
		qprob <- 1 - prob
		lArray <- c(-lambda, -lambda + log(lambda*qprob)) 
		for(x in 2:max(2, xMax)){
			nextProb <- ((lambda*qprob + 2*prob*(x - 1)) - 
							prob^2*(x - 2)*exp(lArray[x - 1] - lArray[x]))/x
			lArray <- c(lArray, log(nextProb) + lArray[x])
			}
		return(lArray)
		}
#
########################################################################################
#
#	Array whose values are the log of the cumulative distribution function, 
#		given an array lArray of values of the log of the probability 
#		function, log(Pr(X=x)).  This is designed to enable one to calculate the 
#		log of the cumulative distribution from the log of the probability when the 
#		probability rounds to zero to machine accuracy.     
#		(note that the array index is 1 more than x) 
#
	gArray <- function(lArray){
		xMax <- length(lArray)
		gArray <- c(lArray[1]) 
		for(x in 1:xMax){
			nextGArray <- gArray[x - 1] + log1p(exp(lArray[x] - gArray[x - 1]))
			gArray <- c(gArray, nextGArray) 
			}
		return(gArray)
		}
#
########################################################################################
#
#	Calculate the log of the tail of the PolyaAeppli distribution, log(Prob(X > x)) 
#		given x, lambda, prob
#
	logTailPA <- function(x, lambda, prob, maxIter=10000){
		xmax <- floor(x)
		last2ells <- dPolyaAeppli(c(xmax, xmax + 1), lambda, prob, log=TRUE)
		lAtXminus2 <- last2ells[1]
		lAtXminus1 <- lAtXmaxPlus1 <- last2ells[2]
		i <- xmax + 2
		sumOfExponentials <- 1
		nextTerm <- 1  # this is to enable while loop to start
		qprob <- 1 - prob
		while(abs(nextTerm) > 2*.Machine$double.eps){
			lAtX <- lAtXminus1 + 
					  log((lambda*qprob + 2*prob*(i - 1) - 
						prob^2*(i - 2)*exp(lAtXminus2 - lAtXminus1))/i)
			nextTerm <- exp(lAtX - lAtXmaxPlus1)
			sumOfExponentials <- sumOfExponentials + nextTerm
#
			lAtXminus2 <- lAtXminus1
			lAtXminus1 <- lAtX
			i <- i + 1
			if(i - xmax - 2 > maxIter){
				warning("\n maxIter exceeded")
				break
				} 
			}
		logTail <- lAtXmaxPlus1 + log(sumOfExponentials)
		return(logTail)
		}
#
########################################################################################
#
#	Array whose values are the log of the upper tail log(Pr(X>x)) of the cumulative 
#		distribution function, given an array lArray of values of the log of the 
#		probability function, log(Pr(X=x)), up to some xmax and the upper tail from 
#		xmax.  This is designed to enable one to 
#		calculate the log of the upper tail of the cumulative distribution from 
#		the log of the probability when the probability rounds to 1 to machine accuracy.     
#		(note that the array index is 1 more than x) 
#
	hArray <- function(hTop, lArray){
		xMax <- length(lArray) - 1
		hArray <- c(hTop)
		for(x in xMax:1){	# calculating hArray[x] = h(x - 1) = log(Pr(X>(x-1)))
			previousHArray <- hArray[1] + log1p(exp(lArray[x + 1] - hArray[1]))
			hArray <- c(previousHArray, hArray) 
			}
		return(hArray)
		}
#
########################################################################################
#
#	Function to detect whole numbers needed for above functions
#
	is.wholenumber <- 
			function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
#
########################################################################################

