# Calculations of effect size estimates of psychedelic-assisted mental health treatments


# General functions
# Make CDF from quantiles
cdfq <- function(q) {
	u <- c(0,.001,.01,.05,.1,.25,.5,.75,.9,.95,.99,.999,1)
	cdf <- spline(q,u, method = 'hyman', xout = seq(-4,4,.001),ties=mean)
	cdf <- unlist(cdf[2], use.names=FALSE)
	return(cdf)
	}

# PDF from quantiles
pdf <- function(q) {
	pdf <- cdfq(q)[1]
	for (i in 2:8001) {
		pdf <- append(pdf,cdfq(q)[i]-cdfq(q)[i-1])
		}
	return(pdf)
	}
	

# Make posterior function, given prior (as a distribution or as quantiles), normal likelihood
post <- function(pr,mu,sig) {
	x <- seq(-4,4,.001)
	lik <- dnorm(x, mean=mu,sd=sig)
	pst <- lik*pr
	pst <- pst/sum(pst)
	return(pst)
	}

# Make discounted posterior via increased variance; faster, if less fine-grained; can reduce search space
postV <- function(pr,mu,sig,d) {
	x <- seq(-4,4,.001)
	m <- sum(pr*x)*d + sum(post(pr,mu,sig)*x)*(1-d)
	t <- c()
	for (i in seq(.001,3,.001)) {
		pst <- post(pr,mu,i)
		t <- append(t,(m - sum(pst*x))^2)
		}
	sd <- which.min(t)*.001
	return(post(pr,mu,sd))
	}
	

# Make discounted posterior via Jeffrey's rule [unused]
postB <- function(pr,mu,sig,pd) {
	u <- c(0,1/4,1/2,3/4,1) # Vector of possible effect sizes, as a fraction of reported effect
	cdfd <- spline(u,pd, method = 'hyman', xout = seq(0,1,.001),ties=mean) # Interpolate to get CDF
	cdfd <- unlist(cdfd[2], use.names=FALSE)
	pdfd <- cdfd[1] # Make pdf
	for (i in 1:1000) {
		pdfd <- append(pdfd,cdfd[i+1] - cdfd[i])
		}
	v <- c() # This will be the matrix of posteriors
	for (i in seq(0,1,.001)) {
		v <- rbind(v,post(pr,i*mu,sig)) # Vector of posteriors
		}
	postB <- pdfd %*% v # Aggregate possible posteriors
	return(postB)
	}

# Calculate moments
meanp <- function(p) {
	x <- seq(-4,4,.001)
	mean <- sum(p*x)
	return(mean)
	}

varp <- function(p) {
	x <- seq(-4,4,.001)
	return(sum(p*(x-meanp(p))^2))
	}

skewp <- function(p) {
	x <- seq(-4,4,.001)
	mu <- meanp(p)
	sig <- sqrt(varp(p))
	return(sum(p*((x-mu)/sig)^3))
	}

kurtp <- function(p) {
	x <- seq(-4,4,.001)
	mu <- meanp(p)
	sig <- sqrt(varp(p))
	return(sum(p*((x-mu)/sig)^4))
	}

# CDF
cdf <- function(p) {
	cdf <- c()
	for (i in 1:8001) {
		cdf <- append(cdf,sum(p[1:i]))
		}
	return(cdf)
	}	
	
# PDF conditional on effect d >= 0.3
pdfcond <- function(p) {
	pc <- p[4301:8001]
	pc <- pc/sum(pc)
	return(pc)
	}

# Conditional mean
meanpc <- function(p) {
	xc <- seq(.3,4,.001)
	return(sum(p*xc))
	}

# CDFs for sensitivity analysis
# Optimistic prior
maxcdfq <- function(q) {
	u <- c(0,.001,.01,.05,.1,.25,.5,.75,.9,.95,.99,.999,1)
	cdf <- rep(u[1],length(seq(q[1]+.001,q[2]-.001,.001)))
	for (i in 2:12) {
		for (j in seq(q[i],q[i+1]-.001,.001)) {
			cdf <- append(cdf,u[i])
			}
		}
	cdf <- append(cdf,c(1,1))
	return(cdf)
	}

# Conservative prior
mincdfq <- function(q) {
	u <- c(0,.001,.01,.05,.1,.25,.5,.75,.9,.95,.99,.999,1)
	cdf <- c(0)
	for (i in 1:12) {
		for (j in seq(q[i]+.001,q[i+1],.001)) {
			cdf <- append(cdf,u[i+1])
			}
		}
	return(cdf)
	}

# PDF from CDF
pdfc <- function(cdf) {
	pdf <- cdf[1]
	for (i in 1:8000) {
		pdf <- append(pdf,cdf[i+1]-cdf[i])
		}
	return(pdf)
	}



# Psilocybin Bayesian Analysis
# Median quantiles
medqP <- c(-4,-1.00,-0.30,-0.10,0.00,0.10,0.18,0.30,0.60,0.80,1.05,1.90,4)

# Prior
priorP <- pdf(medqP)
meanp(priorP) #prior mean
varp(priorP) #prior variance
sum(priorP[4301:8001]) #prob d >= 0.3

# Likelihoods
# Griffiths et al.
muG <- 1.3
sigG <- .381
dG <- .71

# Carhart-Harris et al.
muC <- 1.01
sigC <- .331
dC <- .81

# Posteriors
postG <- postV(priorP,muG,sigG,dG)
postGC <- postV(postG,muC,sigC,dC)

meanp(postG)
sum(postG[4301:8001])

meanp(postGC)
varp(postGC)
sum(postGC[4301:8001])

# Posterior cond. on d >= 0.3
pstcP <- pdfcond(postGC)
meanpc(pstcP)

# 90% CIs
xc <- seq(.3,4,.001)
cdf.pstcP <- cdf(pstcP)
lo <- which.min((cdf.pstcP-.05)^2)
hi <- which.min((cdf.pstcP-.95)^2)
lo <- xc[lo]
hi <- xc[hi]
c(lo,hi)


# Sensitivity analysis
# Conservative
priormin <- mincdfq(medqP)
priormin <- pdfc(priormin)

meanp(priormin) #prior mean
varp(priormin) #prior variance
sum(priormin[4301:8001]) #prob d >= 0.3

postming <- post(priormin,muG,sigG)
postminG <- postV(priormin,muG,sigG,dG)

postminGc <- post(postminG,muC,sigC)
postminGC <- postV(postminG,muC,sigC,dC)

meanp(postming)
sum(postming[4301:8001])
meanp(postminG)
sum(postminG[4301:8001])

meanp(postminGc)
sum(postminGc[4301:8001])

meanp(postminGC)
varp(postminGC)
sum(postminGC[4301:8001])

pstminc <- pdfcond(postminGC)
meanpc(pstminc)

# Optimistic
priormax <- maxcdfq(medqP)
priormax <- pdfc(priormax)

meanp(priormax) #prior mean
varp(priormax) #prior variance
sum(priormax[4301:8001]) #prob d >= 0.3

postmaxg <- post(priormax,muG,sigG)
postmaxG <- postV(priormax,muG,sigG,dG)

postmaxGc <- post(postmaxG,muC,sigC)
postmaxGC <- postV(postmaxG,muC,sigC,dC)

meanp(postmaxg)
sum(postmaxg[4301:8001])
meanp(postmaxG)
sum(postmaxG[4301:8001])

meanp(postmaxGc)
sum(postmaxGc[4301:8001])

meanp(postmaxGC)
varp(postmaxGC)
sum(postmaxGC[4301:8001])

pstmaxc <- pdfcond(postmaxGC)
meanpc(pstmaxc)



# MDMA Bayesian Analysis
# Quantiles
medqM <- c(-4,-0.50,-0.20,-0.05,0.05,0.13,0.25,0.50,0.85,1.15,1.50,1.80,4)

# Prior
priorM <- pdf(medqM)
meanp(priorM) #prior mean
varp(priorM) #prior variance
sum(priorM[4301:8001]) #prob d >= 0.3

# Likelihoods
# Mithoefer et al.
muM1 <- 2.8
sigM1 <- .816
muM2 <- 1.1
sigM2 <- .520
dM <- .74

# Ot'alora et al.
muO1 <- .42
sigO1 <- .529
muO2 <- .37
sigO2 <- .529
dO <- .68

# Posteriors
postM1B <- postV(priorM,muM1,sigM1,dM)
postMB <- postV(postM1B,muM2,sigM2,dM)
postMBO1B <- postV(postMB,muO1,sigO1,dO)
postMBOB <- postV(postMBO1B,muO2,sigO2,dO)

meanp(postMB)
sum(postMB[4301:8001])

meanp(postMBOB)
varp(postMBOB)
sum(postMBOB[4301:8001])

#posterior cond. on d >= 0.3
pstcM <- postMBOB[4301:8001]
pstcM <- pstcM/sum(pstcM)
sum(pstcM*seq(.3,4,.001))

# 90% CIs
xc <- seq(.3,4,.001)
cdf.pstcM <- cdf(pstcM)
lo <- which.min((cdf.pstcM-.05)^2)
hi <- which.min((cdf.pstcM-.95)^2)
lo <- xc[lo]
hi <- xc[hi]
c(lo,hi)

# Sensitivity analysis
# Conservative
priormin <- mincdfq(medqM)
priormin <- pdfc(priormin)

meanp(priormin) #prior mean
varp(priormin) #prior variance
sum(priormin[4301:8001]) #prob d >= 0.3

postminMV <- postV(postV(priormin,muM1,sigM1,dM),muM2,sigM2,dM)
postminMVOV <- postV(postV(postminMV,muO1,sigO1,dO),muO2,sigO2,dO)

meanp(postminMV)
sum(postminMV[4301:8001])

meanp(postminMVOV)
varp(postminMVOV)
sum(postminMVOV[4301:8001])

pstminc <- postminMVOV[4301:8001]
pstminc <- pstminc/sum(pstminc)
sum(pstminc*seq(.3,4,.001))

# Optimistic
priormax <- maxcdfq(medqM)
priormax <- pdfc(priormax)

meanp(priormax) #prior mean
varp(priormax) #prior variance
sum(priormax[4301:8001]) #prob d >= 0.3

postmaxMV <- postV(postV(priormax,muM1,sigM1,dM),muM2,sigM2,dM)
postmaxMVOV <- postV(postV(postmaxMV,muO1,sigO1,dO),muO2,sigO2,dO)

meanp(postmaxMV)
sum(postmaxMV[4301:8001])

meanp(postmaxMVOV)
varp(postmaxMVOV)
sum(postmaxMVOV[4301:8001])

pstmaxc <- postmaxMVOV[4301:8001]
pstmaxc <- pstmaxc/sum(pstmaxc)
sum(pstmaxc*seq(.3,4,.001))




# Non-Bayesian effect size estimates

# PDF from 5 points on CDF, as a proportion of reported effect size
prob <- function(p) {
	u <- c(0,1/4,1/2,3/4,1) # Vector of possible effect sizes, as a fraction of reported effect
	cdfd <- spline(u,p, method = 'hyman', xout = seq(0,1,.001),ties=mean) # Interpolate to get CDF
	cdfd <- unlist(cdfd[2], use.names=FALSE)
	pdfd <- cdfd[1] # Make pdf
	for (i in 1:1000) {
		pdfd <- append(pdfd,cdfd[i+1] - cdfd[i])
		}
	return(pdfd)
	}

prob2 <- function(p) {
	u <- c(-1,0,1/4,1/2,3/4,1) # Vector of possible effect sizes, as a fraction of reported effect
	cdfd <- spline(u,p, method = 'hyman', xout = seq(-1,1,.001),ties=mean) # Interpolate to get CDF
	cdfd <- unlist(cdfd[2], use.names=FALSE)
	pdfd <- cdfd[1] # Make pdf
	for (i in 1:2000) {
		pdfd <- append(pdfd,cdfd[i+1] - cdfd[i])
		}
	return(pdfd)
	}



# Psilocybin
# Reported effect size
muP <- .83

# 5 points on effect size CDF, as a proportion of reported effect size
pP <- c(0,.55,.865,.965,1)

# Effect size PDF
probP <- prob(pP)

# x-coordinates between 0 and reported effect size
xP <- seq(0,1,.001)
xP <- muP*xP

# Results
meanP <- sum(xP*probP)
meanP
sum(probP[363:1001])

probPc <- probP[363:1001]
probPc <- probPc/sum(probPc)
sum(probPc*xP[363:1001])

cdf.probPc <- cdf(probPc)
lo <- which.min((cdf.probPc-.05)^2)
hi <- which.min((cdf.probPc-.95)^2)
lo <- xP[363:1001][lo]
hi <- xP[363:1001][hi]
c(lo,hi)





# MDMA
# Reported effect size
muM <- .8

# 5 points on effect size CDF, as a proportion of reported effect size
pM <- c(0,.515,.83,.95,1)

# Effect size PDF
probM <- prob(pM)

# x-coordinates between 0 and reported effect size
xM <- seq(0,1,.001)
xM <- muM*xM

# Results
meanM <- sum(xM*probM)
meanM
sum(probM[376:1001])

probMc <- probM[376:1001]
probMc <- probMc/sum(probMc)
sum(probMc*xM[376:1001])

cdf.probMc <- cdf(probMc)
lo <- which.min((cdf.probPc-.05)^2)
hi <- which.min((cdf.probPc-.95)^2)
lo <- xM[376:1001][lo]
hi <- xM[376:1001][hi]
c(lo,hi)



# Graphs
x <- seq(-4,4,.001)

# Psilocybin
# Prior
plot(x[3500:5500],priorP[3500:5500], type="l", col="blue", lwd=1, xlab="Effect size", ylab="Probability density", main="Psilocybin prior PDF")

# Posterior
plot(x[3500:6000],postGC[3500:6000], type="l", col="red", lwd=1, xlab="Effect size", ylab="Probability density", main="Psilocybin posterior PDF")

# Prior & posterior
plot(x[3500:6000],priorP[3500:6000], type="l", col="blue", lwd=1, xlab="Effect size", ylab="Probability density", main="Psilocybin prior and posterior PDFs")
lines(x[3500:6000],postGC[3500:6000], type="l", col="red", lwd=1)
legend(x = .6, y = .003, c("Psilocybin prior","Psilocybin posterior"), lty = c(1,1), lwd = c(2,2), col = c("blue","red"), bg = 'white')

# Non-Bayesian PDF
plot(xP,probP, type="l", col="orange", lwd=1, xlab="Effect size", ylab="Probability density", main="Psilocybin Non-Bayesian PDF")



# MDMA
# Prior
plot(x[3500:6000],priorM[3500:6000], type="l", col="turquoise", lwd=1, xlab="Effect size", ylab="Probability density", main="MDMA prior PDF")

# Posterior
plot(x[3500:6000],postMBOB[3500:6000], type="l", col="dark red", lwd=1, xlab="Effect size", ylab="Probability density", main="MDMA posterior PDF")

# Prior & posterior
plot(x[3500:6000],priorM[3500:6000], type="l", col="turquoise", lwd=1, xlab="Effect size", ylab="Probability density", main="MDMA prior and posterior PDFs")
lines(x[3500:6000],postMBOB[3500:6000], type="l", col="dark red", lwd=1)
legend(x = .6, y = .002, c("MDMA prior","MDMA posterior"), lty = c(1,1), lwd = c(2,2), col = c("turquoise","dark red"), bg = 'white')

# Non-Bayesian PDF
plot(xM,probM, type="l", col="brown", lwd=1, xlab="Effect size", ylab="Probability density", main="MDMA Non-Bayesian PDF")
