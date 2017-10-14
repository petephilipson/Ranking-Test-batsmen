# File created 09/08/17 for quadratic ageing (LN prior for curvature)
batFit <- function(data, niter, thin, burnin, mudecpriors = NULL, sigpriors = NULL, sigdelpriors = NULL, alphapriors = NULL, hetpriors = NULL, pizippriors = NULL, corner1 = FALSE, inits = NULL, assoc_del = 1)
{
library(MASS)
mc <- match.call()
# Get years and age indexed from 1 upwards
YearScale <- data$Year - (min(data$Year) - 1)
AgeScale <- data$Age - (min(data$Age) - 1) # Mean correct this?
data <- cbind(data, YearScale, AgeScale)
# Append with some factors which save time when using tapply
data$IDF <- as.factor(data$ID)
data$YearScaleF <- as.factor(data$YearScale)
data$AgeScaleF <- as.factor(data$AgeScale)
# Index for not-outs
censInd <- (data$NO == 1)
# Index for zero-inflation
zeroInd <- (data$Runs == 0)
D <- sum(zeroInd)
# Get index data
np <- length(unique(data$ID)) # No. of players
ny <- max(YearScale) - min(YearScale) + 1 # No. of years (includes wars - might remove later)
na <- max(AgeScale) - min(AgeScale) + 1 # No. of age parameters
nd <- max(data$Decade) # No. of birth decades 
nera <- max(data$Era)
nopp <- length(table(data$Opposition))
ninns <- 4
# Set up result matrix
npara <- 2*nd + 5*np + ny + ninns + nopp*nera + 5
matres <- matrix(0, ncol = 1 + npara, nrow = niter)
matres[, 1] <- 1:niter
# Set up fitted values matrix
fitres <- matrix(0, ncol = 2*dim(data)[1] + 1, nrow = niter) 
fitres[, 1] <- 1:niter
# Define penalty matrix for random walk block updates for year effects
# Year effects: first order random walk
penM <- diag(ny)
diag(penM) <- diag(penM) * c(1, rep(1 + assoc_del^2, ny - 2), 1)
diag(penM[1:(ny - 1), 2:ny]) <- -assoc_del
diag(penM[2:ny, 1:(ny - 1)]) <- -assoc_del
if (corner1) {penM <- penM[-1, -1]}
else {penM <- penM[-ny, -ny]}
# Get summary data (from function 'summarydata')
summ <- summarydata(data)
if(is.null(mudecpriors)){mudecpriors <- c(log(20), 0.25)} 
if(is.null(sigpriors)){sigpriors <- c(3, 1)} 
if(is.null(sigdelpriors)){sigdelpriors <- c(2, 0.01)} 
if(is.null(alphapriors)){alphapriors <- c(30, 2, -6, 2)}
if(is.null(hetpriors)){hetpriors <- c(0, 0.5)}
if(is.null(pizippriors)){pizippriors <- c(0.5, 4.5)}
paras <- initialise(data, inits, summ$yInd, np, ny, nd, nera, na, ninns, nopp, pdec = summ$pdec, k = summ$k, penM, mudecpriors, sigpriors, sigdelpriors, alphapriors, hetpriors, pizippriors)
res <- NULL
for (i in 1 : burnin){
	print(i)
	poisparas <- lambdaUpdate(data, nera, censInd, zeroInd, paras$eta, paras$lambda, paras$theta, paras$delta, paras$agefun, paras$xi, paras$kappa, paras$omega, paras$phi, paras$pizip, rwsd = 1)
	decparas <- mudecsigUpdate(nd, summ$pdec, summ$k, mudecpriors, sigpriors, mudec = paras$mudec, sig = paras$sig, thbar = paras$thetabar, paras$theta)
	playerparas <- thetaUpdate(data, np, nera, summ$pdec, summ$pinns, decparas$mup, decparas$sigp, paras$eta, poisparas$lambda, paras$theta, paras$delta, paras$agefun, paras$xi, paras$kappa, paras$omega, paras$phi, rwsd = 0.3)
	yearparas <- deltaUpdate(data, ny, nera, summ$yyears, summ$yInd, penM, sigdelpriors, paras$sigdel, paras$eta, poisparas$lambda, playerparas$theta, paras$delta, paras$agefun, paras$xi, paras$kappa, paras$omega, paras$phi, corner1)
	ageparas <- alphaUpdate(data, np, nera, paras$mual1, paras$mual2, paras$sigal1, paras$sigal2, paras$eta, poisparas$lambda, playerparas$theta, yearparas$delta, paras$agefun, paras$alphamat, paras$xi, paras$kappa, paras$omega, paras$phi, rwsd1 = 0.85, rwsd2 = 3)
	innsparas <- innsUpdate(data, nera, paras$sigi, paras$eta, poisparas$lambda, playerparas$theta, yearparas$delta, ageparas$agefun, paras$xi, paras$kappa, paras$omega, paras$phi, rwsd = 0.04)
	awaypara <- awayUpdate(data, nera, paras$sigk, paras$eta, poisparas$lambda, playerparas$theta, yearparas$delta, ageparas$agefun, innsparas$xi, paras$kappa, paras$omega, paras$phi, rwsd = 0.02)
	handpara <- handUpdate(data, nera, paras$sigo, paras$eta, poisparas$lambda, playerparas$theta, yearparas$delta, ageparas$agefun, innsparas$xi, awaypara$kappa, paras$omega, paras$phi, rwsd = 0.04)
	oppparas <- oppUpdate(data, nopp, nera, summ$oppdataInd, paras$sigopp, paras$eta, poisparas$lambda, playerparas$theta, yearparas$delta, ageparas$agefun, innsparas$xi, awaypara$kappa, handpara$omega, paras$phi, rwsd = 0.35)
	hetpara <- etaUpdate(data, np, nera, summ$pinns, hetpriors, paras$eta, poisparas$lambda, playerparas$theta, yearparas$delta, ageparas$agefun, innsparas$xi, awaypara$kappa, handpara$omega, oppparas$phi, rwsd = 1)
	pizippara <- pizipUpdate(data, np, zeroInd, D, censInd, summ$pinns, summ$pducks, paras$pizipa, paras$pizipb, poisparas$lambda, paras$pizip, rwsd = 2)
	paras <- list(mudec = decparas$mudec, sig = decparas$sig, mup = decparas$mup, sigp = decparas$sigp, theta = playerparas$theta, thetabar = playerparas$thetabar, lambda = poisparas$lambda, eta = hetpara$eta, etal = hetpara$etal, delta = yearparas$delta, sigdel = yearparas$sigdel, alphamat = ageparas$alphamat, alpha1 = ageparas$alpha1, alpha2 = ageparas$alpha2, mual1 = ageparas$mual1, mual2 = ageparas$mual2, sigal1 = ageparas$sigal1, sigal2 = ageparas$sigal2, agefun = ageparas$agefun, xi = innsparas$xi, sigi = innsparas$sigi, kappa = awaypara$kappa, sigk = awaypara$sigk, omega = handpara$omega, sigo = handpara$sigo, phi = oppparas$phi, sigopp = oppparas$sigopp, pizip = pizippara$pizip, pizipa = pizippara$pizipa, pizipb = pizippara$pizipb, mu = hetpara$mu)
}
for (j in 1 : niter){
	print(j)
	for (k in 1 : thin){
		poisparas <- lambdaUpdate(data, nera, censInd, zeroInd, paras$eta, paras$lambda, paras$theta, paras$delta, paras$agefun, paras$xi, paras$kappa, paras$omega, paras$phi, paras$pizip, rwsd = 1)
		decparas <- mudecsigUpdate(nd, summ$pdec, summ$k, mudecpriors, sigpriors, mudec = paras$mudec, sig = paras$sig, thbar = paras$thetabar, paras$theta)
		playerparas <- thetaUpdate(data, np, nera, summ$pdec, summ$pinns, decparas$mup, decparas$sigp, paras$eta, poisparas$lambda, paras$theta, paras$delta, paras$agefun, paras$xi, paras$kappa, paras$omega, paras$phi, rwsd = 0.3)
		yearparas <- deltaUpdate(data, ny, nera, summ$yyears, summ$yInd, penM, sigdelpriors, paras$sigdel, paras$eta, poisparas$lambda, playerparas$theta, paras$delta, paras$agefun, paras$xi, paras$kappa, paras$omega, paras$phi, corner1)
	ageparas <- alphaUpdate(data, np, nera, paras$mual1, paras$mual2, paras$sigal1, paras$sigal2, paras$eta, poisparas$lambda, playerparas$theta, yearparas$delta, paras$agefun, paras$alphamat, paras$xi, paras$kappa, paras$omega, paras$phi, rwsd1 = 0.85, rwsd2 = 3)
		innsparas <- innsUpdate(data, nera, paras$sigi, paras$eta, poisparas$lambda, playerparas$theta, yearparas$delta, ageparas$agefun, paras$xi, paras$kappa, paras$omega, paras$phi, rwsd = 0.04)
		awaypara <- awayUpdate(data, nera, paras$sigk, paras$eta, poisparas$lambda, playerparas$theta, yearparas$delta, ageparas$agefun, innsparas$xi, paras$kappa, paras$omega, paras$phi, rwsd = 0.02)
		handpara <- handUpdate(data, nera, paras$sigo, paras$eta, poisparas$lambda, playerparas$theta, yearparas$delta, ageparas$agefun, innsparas$xi, awaypara$kappa, paras$omega, paras$phi, rwsd = 0.04)
		oppparas <- oppUpdate(data, nopp, nera, summ$oppdataInd, paras$sigopp, paras$eta, poisparas$lambda, playerparas$theta, yearparas$delta, ageparas$agefun, innsparas$xi, awaypara$kappa, handpara$omega, paras$phi, rwsd = 0.35)
		hetpara <- etaUpdate(data, np, nera, summ$pinns, hetpriors, paras$eta, poisparas$lambda, playerparas$theta, yearparas$delta, ageparas$agefun, innsparas$xi, awaypara$kappa, handpara$omega, oppparas$phi, rwsd = 1)
		pizippara <- pizipUpdate(data, np, zeroInd, D, censInd, summ$pinns, summ$pducks, paras$pizipa, paras$pizipb, poisparas$lambda, paras$pizip, rwsd = 2)
	paras <- list(mudec = decparas$mudec, sig = decparas$sig, mup = decparas$mup, sigp = decparas$sigp, theta = playerparas$theta, thetabar = playerparas$thetabar, lambda = poisparas$lambda, eta = hetpara$eta, etal = hetpara$etal, delta = yearparas$delta, sigdel = yearparas$sigdel, alphamat = ageparas$alphamat, alpha1 = ageparas$alpha1, alpha2 = ageparas$alpha2, mual1 = ageparas$mual1, mual2 = ageparas$mual2, sigal1 = ageparas$sigal1, sigal2 = ageparas$sigal2, agefun = ageparas$agefun, xi = innsparas$xi, sigi = innsparas$sigi, kappa = awaypara$kappa, sigk = awaypara$sigk, omega = handpara$omega, sigo = handpara$sigo, phi = oppparas$phi, sigopp = oppparas$sigopp, pizip = pizippara$pizip, pizipa = pizippara$pizipa, pizipb = pizippara$pizipb, mu = hetpara$mu)
#		loglikelihood <- lhood(data, nera, paras$theta, paras$delta, paras$alpha1, paras$alpha2, paras$alpha3, paras$xi, paras$kappa, paras$omega, paras$phi, paras$lambda, censInd, zeroInd)
	}
	matres[j, -1] <- c(paras$alpha1, paras$alpha2, paras$theta, paras$eta, paras$phi, paras$mudec, paras$sig, paras$delta, paras$xi, paras$kappa, paras$omega, paras$pizip, paras$sigdel)
	fitres[j, -1] <- c(paras$mu, paras$lambda)
}
list(res = matres, fitvals = fitres, dataset = data, rwsds = batFit, inits = initialise, call = mc)
}

# Define age function
agecurve <- function(Age, alphamatl){
agecurve <- -alphamatl[, 2]*(alphamatl[, 1] - Age)^2
agecurve
}

summarydata <- function(data)
{
pinns <- as.numeric(table(data$ID)) # No. of innings played for each player
pruns <- with(data, tapply(Runs, ID, sum)) # Runs per player
pyears <- with(data, tapply(Year, ID, unique)) # Years per player
pducks <- with(data, tapply(Runs == 0, ID, sum))# Ducks per player
aruns <- with(data, tapply(Runs, AgeScale, sum)) # Runs per age
ainns <- with(data, tapply(Inns, AgeScale, length)) # Innings per age
anos <- with(data, tapply(NO, AgeScale, sum)) # Not-outs per age
ydata1 <- as.numeric(table(data$YearScale)) # Data per year - careful as some years have no data
ydata2 <- with(data, tapply(Runs, YearScale, sum)) # Runs per year - careful as some years have no data
ydataInd <- with(data, 1:max(YearScale) %in% unique(YearScale))
yyears <- with(data, rep(0, max(YearScale)))
yruns <- yyears
yyears[ydataInd] <- ydata1
yruns[ydataInd] <- ydata2
decCount <- with(data, tapply(Decade, ID, unique)) # Players born in each decade
pdec <- with(data, tapply(Decade, ID, sum))/pinns # Vector of decades born, one entry for each player
k <- as.numeric(table(unlist(decCount))) # Use k to match previous work, as.numeric is nicer to work with
ave <- pruns/pinns
ageave <- aruns/(ainns - anos)
oppdataInd <- as.vector(!is.na(with(data, tapply(Runs, OppEra), length)))
print(oppdataInd)
list(ave = ave, pyears = pyears, pruns = pruns, pinns = pinns, pducks = pducks, yruns = yruns, yyears = yyears, yInd = ydataInd, pdec = pdec, k = k, oppdataInd = oppdataInd)
}

initialise <- function(data, inits, yInd, np, ny, nd, nera, na, ninns, nopp, pdec, k, penM, mudecpriors, sigpriors, sigdelpriors, alphapriors, hetpriors, pizippriors)
{
# These are all variances	
mual1 <- alphapriors[1]
mual2 <- alphapriors[3]
sigal1 <- alphapriors[2]^2
sigal2 <- alphapriors[4]^2
sigi <- 0.5^2
sigk <- 0.5^2
sigo <- 0.5^2
sigopp <- 0.25^2
pizipa <- pizippriors[1]
pizipb <- pizippriors[2]
# Set up ageing prior
alpha1 <- rnorm(np, mual1, sqrt(sigal1))
alpha2 <- rlnorm(np, mual2, sqrt(sigal2))
# Initialise using supplied starting values or sampling from priors
if(is.null(inits)){
	mudec <- rnorm(nd, mudecpriors[1], mudecpriors[2])
	sigdec <- 1/rgamma(nd, shape = sigpriors[1], rate = sigpriors[2])
	theta <- rnorm(np, mudec[pdec], sigdec[pdec])
	delta <- rep(0, ny)
	xi <- rnorm(ninns, 0, sqrt(sigi))
	kappa <- rnorm(2, 0, sqrt(sigk))
	phi <- rep(rnorm(nopp, 0, sqrt(sigopp)), each = nera)
	omega <- rnorm(2, 0, sqrt(sigo))
	eta <- rlnorm(np, hetpriors[1], hetpriors[2])
	pizip <- rbeta(np, pizipa, pizipb)
}
else {
	theta <- inits$theta0
	mudec <- mean(theta)
	sigdec <- var(theta)
	delta <- rep(0, ny)
	delta[yInd] <- inits$delta0
	xi <- inits$xi0
	kappa <- inits$kappa0 
	omega <- rep(0, 2)
	phi <- rep(inits$phi0, each = nera)
	eta <- rep(inits$eta0, np)
	pizip <- rep(inits$pi0, np)
	}
# Also need thetabar for each decade for mudec and sigma conditionals
thetabar <- tapply(theta, pdec, mean)
# Random walk variance from inverse gamma for years
sigdel <- 1/rgamma(1, shape = sigdelpriors[1], rate = sigdelpriors[2])
# Initialise Poisson parameter
alphamat <- cbind(alpha1, alpha2)
agefun <- with(data, agecurve(Age, alphamat[ID, ]))
mu <- with(data, exp(theta[ID] + delta[YearScale] + agefun  + xi[Inns] + kappa[HA] + omega[Hand] + phi[(Opposition - 1)*nera + Era]))
lambda <- rgamma(length(mu), shape = eta[data$ID], rate = eta[data$ID]/mu)
etal <- eta[data$ID]
list(mudec = mudec, sig = sigdec, theta = theta, thetabar = thetabar, delta = delta, sigdel = sigdel, lambda = lambda, eta = eta, etal = etal, agefun = agefun, alpha1 = alpha1, alpha2 = alpha2, mual1 = mual1, mual2 = mual2, sigal1 = sigal1, sigal2 = sigal2, alphamat = alphamat, xi = xi, sigi = sigi, kappa = kappa, sigk = sigk, phi = phi, sigopp = sigopp, omega = omega, sigo = sigo, pizip = pizip, pizipa = pizipa, pizipb = pizipb)
}

# Heterogeneity parameter update
etaUpdate <- function(data, np, nera, pinns, hetpriors, eta, lambda, theta, delta, agefun, xi, kappa, omega, phi, rwsd){
mueta <- hetpriors[1]
sigeta <- hetpriors[2]
mu <- with(data, exp(theta[ID] + delta[YearScale] + agefun + xi[Inns] + kappa[HA] + omega[Hand] + phi[(Opposition - 1)*nera + Era]))
Ni <- pinns
eta1 <- rlnorm(np, log(eta), rwsd)
denrat <- (mueta/sigeta^2 - 1)*log(eta1/eta) + (log(eta)^2 - log(eta1)^2)/(2*sigeta^2) + Ni*(lgamma(eta) - lgamma(eta1) + eta1*log(eta1) - eta*log(eta))
denrat <- denrat + (eta1 - eta)*with(data, tapply(log(lambda/mu) - lambda/mu, ID, sum)) + log(eta1/eta)
laccept <- pmin(0, denrat)
accept <- (log(runif(np)) < laccept)
eta[accept] <- eta1[accept]
etal <- eta[data$ID]
list(eta = eta, etal = etal, mu = mu)	
}

# Poisson RE parameter update
lambdaUpdate <- function(data, nera, censInd, zeroInd, eta, lambda, theta, delta, agefun, xi, kappa, omega, phi, pizip, rwsd){
mu <- with(data, exp(theta[ID] + delta[YearScale] + agefun + xi[Inns] + kappa[HA] + omega[Hand] + phi[(Opposition - 1)*nera + Era]))
# Combined update
lambda1 <- rlnorm(dim(data)[1], log(lambda), rwsd)
cenprob1 <- with(data, ppois(Runs - 1, lambda1, lower.tail = FALSE, log.p = TRUE))
cenprob <- with(data, ppois(Runs - 1, lambda, lower.tail = FALSE, log.p = TRUE))
denrat <- (eta[data$ID] - 1 + data$Runs*(1 - censInd)*(1 - zeroInd))*log(lambda1/lambda) - (eta[data$ID]/mu + (1 - censInd)*(1 - zeroInd))*(lambda1 - lambda) + censInd*(1 - zeroInd)*(cenprob1 - cenprob) 
denrat <- denrat + log(lambda1/lambda)
denrat <- denrat + with(data, zeroInd*(log(pizip[ID] + (1 - pizip[ID])*exp(-lambda1*(1 - censInd))) - log(pizip[ID] + (1 - pizip[ID])*exp(-lambda*(1 - censInd)))))
laccept <- pmin(0, denrat)
accept <- (log(runif(dim(data)[1])) < laccept)
lambda[accept] <- lambda1[accept]
list(lambda = lambda)
}

# Decade parameter updates
mudecsigUpdate <- function(nd, pdec, k, mudecpriors, sigpriors, mudec, sig, thbar, theta){
m <- mudecpriors[1]
s <- mudecpriors[2]
a <- sigpriors[1]
b <- sigpriors[2]
mean <- (m/(s^2) + k*thbar/sig)/((1/(s^2) + k/sig))
sd <- sqrt(1/(1/(s^2) + k/sig))
mudec <- rnorm(nd, mean, sd)
a <- a + 0.5*k
b <- b + 0.5*tapply((theta - mudec[pdec])^2, pdec, sum)
sig <- 1/rgamma(nd, shape = a, rate = b)
mup <- mudec[pdec]
sigp <- sig[pdec]
list(mudec = mudec, sig = sig, mup = mup, sigp = sigp)
}

# Player parameter updates
thetaUpdate <- function(data, np, nera, pdec, pinns, mup, sigp, eta, lambda, theta, delta, agefun, xi, kappa, omega, phi, rwsd){
theta1 <- rnorm(np, theta, rwsd)
denrat <- ((theta - mup)^2 - (theta1 - mup)^2)/(2*sigp) - eta*(exp(-theta1) - exp(-theta))*with(data, tapply(lambda*exp(-(delta[YearScale] + agefun + xi[Inns] + kappa[HA] + omega[Hand] + phi[(Opposition - 1)*nera + Era])), ID, sum))
denrat <- denrat - eta*pinns*(theta1 - theta)
laccept <- pmin(0, denrat)
accept <- (log(runif(np)) < laccept)
theta[accept] <- theta1[accept]
#print(summary(theta))
thetabar <- tapply(theta, pdec, mean)
#cat("theta:", round(100*sum(accept)/np, 0), "% ")
list(theta = theta, thetabar = thetabar)
}  

# Year parameter updates
deltaUpdate <- function(data, ny, nera, yyears, yind, penM, sigdelpriors, sigdel, eta, lambda, theta, delta, agefun, xi, kappa, omega, phi, corner1){
if (corner1){
delta <- delta[-1]
yind <- yind[-1]
remove <- data$YearScale == 1
}
else{
delta <- delta[-ny]
yind <- yind[-ny]
remove <- data$YearScale == max(data$YearScale)
}
# Remove data from corner-constrained year
dataC <- data[!remove, ]
dataC$YearScale <- dataC$YearScale - min(dataC$YearScale) + 1
lambdaC <- lambda[!remove]
agefunC <- agefun[!remove]
# Careful with years without data
Q <- penM/sigdel
# Now proposal distribution - needs approximations
# First approximation
b0 <- rep(0, ny - 1)
c0 <- rep(0, ny - 1)
# Careful - years without data
b0[yind] <- with(dataC, tapply(eta[ID]*(lambdaC*exp(- theta[ID] - delta[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era])*(1 + delta[YearScale]) - 1), YearScale, sum))
c0[yind] <- with(dataC, tapply(eta[ID]*lambdaC*exp(- theta[ID] - delta[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era]), YearScale, sum))
Q0 <- Q + diag(c0)
mu1 <- solve(Q0, b0)
# delta1 <- mvrnorm(n = 1, mu = mu1, Sigma = solve(Q0)) 
# Second approximation
b1 <- rep(0, ny - 1)
c1 <- rep(0, ny - 1)
b1[yind] <- with(dataC, tapply(eta[ID]*(lambdaC*exp(- theta[ID] - mu1[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era])*(1 + mu1[YearScale]) - 1), YearScale, sum))
c1[yind] <- with(dataC, tapply(eta[ID]*lambdaC*exp(- theta[ID] - mu1[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era]), YearScale, sum))
Q1 <- Q + diag(c1)
mu2 <- solve(Q1, b1)
# delta2 <- mvrnorm(n = 1, mu = mu2, Sigma = solve(Q1)) 
# Third approximation
b2 <- rep(0, ny - 1)
c2 <- rep(0, ny - 1)
b2[yind] <- with(dataC, tapply(eta[ID]*(lambdaC*exp(- theta[ID] - mu2[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era])*(1 + mu2[YearScale]) - 1), YearScale, sum))
c2[yind] <- with(dataC, tapply(eta[ID]*lambdaC*exp(- theta[ID] - mu2[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era]), YearScale, sum))
Q2 <- Q + diag(c2)
mu3 <- solve(Q2, b2)
# delta3 <- mvrnorm(n = 1, mu = mu3, Sigma = solve(Q2))
# Fourth approximation
b3 <- rep(0, ny - 1)
c3 <- rep(0, ny - 1)
b3[yind] <- with(dataC, tapply(eta[ID]*(lambdaC*exp(- theta[ID] - mu3[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era])*(1 + mu3[YearScale]) - 1), YearScale, sum))
c3[yind] <- with(dataC, tapply(eta[ID]*lambdaC*exp(- theta[ID] - mu3[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era]), YearScale, sum))
Q3 <- Q + diag(c3)
mu4 <- solve(Q3, b3)
delta4 <- mvrnorm(n = 1, mu = mu4, Sigma = solve(Q3))
# Define final term
deltanew <- delta4
mucurr <- mu4
Qcurr <- Q3
# Construct proposal ratio term
# First approximation for q(del | delstar)
bnew <- rep(0, ny - 1)
cnew <- rep(0, ny - 1)
bnew[yind] <- with(dataC, tapply(eta[ID]*(lambdaC*exp(- theta[ID] - deltanew[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era])*(1 + deltanew[YearScale]) - 1), YearScale, sum))
cnew[yind] <- with(dataC, tapply(eta[ID]*lambdaC*exp(- theta[ID] - deltanew[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era]), YearScale, sum))
Qnew <- Q + diag(cnew)
munew1 <- solve(Qnew, bnew)
# Second approximation for q(del | delstar)
bnew1 <- rep(0, ny - 1)
cnew1 <- rep(0, ny - 1)
bnew1[yind] <- with(dataC, tapply(eta[ID]*(lambdaC*exp(- theta[ID] - munew1[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era])*(1 + munew1[YearScale]) - 1), YearScale, sum))
cnew1[yind] <- with(dataC, tapply(eta[ID]*lambdaC*exp(- theta[ID] - munew1[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era]), YearScale, sum))
Qnew1 <- Q + diag(cnew1)
munew2 <- solve(Qnew1, bnew1)
# Third approximation for q(del | delstar)
bnew2 <- rep(0, ny - 1)
cnew2 <- rep(0, ny - 1)
bnew2[yind] <- with(dataC, tapply(eta[ID]*(lambdaC*exp(- theta[ID] - munew2[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era])*(1 + munew2[YearScale]) - 1), YearScale, sum))
cnew2[yind] <- with(dataC, tapply(eta[ID]*lambdaC*exp(- theta[ID] - munew2[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era]), YearScale, sum))
Qnew2 <- Q + diag(cnew2)
munew3 <- solve(Qnew2, bnew2)
# Fourth approximation for q(del | delstar)
bnew3 <- rep(0, ny - 1)
cnew3 <- rep(0, ny - 1)
bnew3[yind] <- with(dataC, tapply(eta[ID]*(lambdaC*exp(- theta[ID] - munew3[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era])*(1 + munew3[YearScale]) - 1), YearScale, sum))
cnew3[yind] <- with(dataC, tapply(eta[ID]*lambdaC*exp(- theta[ID] - munew3[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era]), YearScale, sum))
Qnew3 <- Q + diag(cnew3)
munew4 <- solve(Qnew3, bnew3)
# Construct proposal ratio
muprop <- munew4
Qprop <- Qnew3
proprat <- 0.5*(log(det(solve(Qcurr, Qprop))) + (deltanew - mucurr)%*%Qcurr%*%(deltanew - mucurr) - (delta - muprop)%*%Qprop%*%(delta - muprop))
# Now density ratio using proposed values 'deltanew'
denrat <- 0.5*(delta%*%Q%*%delta - deltanew%*%Q%*%deltanew) - with(dataC, sum(eta[ID]*(lambdaC*(exp(- theta[ID] - deltanew[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era]) - exp(- theta[ID] - delta[YearScale] - agefunC - xi[Inns] - kappa[HA] - omega[Hand] - phi[(Opposition - 1)*nera + Era])) + deltanew[YearScale] - delta[YearScale])))
laccept <- pmin(0, denrat + proprat)
accept <- (log(runif(1)) < laccept)
if (accept) {delta <- deltanew}
# Update variance term
g <- sigdelpriors[1] + 0.5*(ny - 1)
h <- sigdelpriors[2] + 0.5*(delta%*%penM%*%delta)
sigdel <- 1/rgamma(1, shape = g, rate = h)
if (corner1){delta <- c(0, delta)}
else {delta <- c(delta, 0)}
# delta <- rep(0, ny)
# print(summary(delta))
list(delta = delta, sigdel = sigdel)
}

# Age parameter updates for quadratic individual ageing with log-normal prior for alpha2
alphaUpdate <- function(data, np, nera, mual1, mual2, sigal1, sigal2, eta, lambda, theta, delta, agefun, alphamat, xi, kappa, omega, phi, rwsd1, rwsd2){
alpha1 <- alphamat[, 1]
alpha2 <- alphamat[, 2]
alpha1c <- rnorm(np, alpha1, rwsd1)
alpha2c <- rlnorm(np, log(alpha2), rwsd2)
alphamat_c <- cbind(alpha1c, alpha2c)
# Define density ratios, acceptance etc.
agefunc <- with(data, agecurve(Age, alphamat_c[ID, ]))
denrat <- ((alpha1 - mual1)^2 - (alpha1c - mual1)^2)/(2*sigal1) + (mual2/sigal2 - 1)*log(alpha2c/alpha2) + (log(alpha2)^2 - log(alpha2c)^2)/(2*sigal2) - eta*(exp(-theta))*with(data, tapply(lambda*exp(-(delta[YearScale] + xi[Inns] + kappa[HA] + omega[Hand] + phi[(Opposition - 1)*nera + Era]))*(exp(-agefunc) - exp(-agefun)), ID, sum))
denrat <- denrat - with(data, tapply(eta[ID]*(agefunc - agefun), ID, sum)) + log(alpha2c/alpha2)
laccept <- pmin(0, denrat)
accept <- (log(runif(np)) < laccept)
alphamat[accept, ] <- alphamat_c[accept, ]
agefun <- with(data, agecurve(Age, alphamat[ID, ]))
list(agefun = agefun, alphamat = alphamat, alpha1 = alphamat[, 1], alpha2 = alphamat[, 2], mual1 = mual1, mual2 = mual2, sigal1 = sigal1, sigal2 = sigal2)
}

# Innings effect parameter update(s)
innsUpdate <- function(data, nera, sigi, eta, lambda, theta, delta, agefun, xi, kappa, omega, phi, rwsd){
xi1 <- rnorm(4, xi, rwsd)
denrat <- (xi^2 - xi1^2)/(2*sigi) - (exp(-xi1) - exp(-xi))*with(data, tapply(eta[ID]*lambda*exp(-(theta[ID] + delta[YearScale] + agefun + kappa[HA] + omega[Hand] + phi[(Opposition - 1)*nera + Era])), Inns, sum))
denrat <- denrat - (xi1 - xi)*with(data, tapply(eta[ID], Inns, sum))
laccept <- pmin(0, denrat)
accept <- (log(runif(4)) < laccept)
xi[accept] <- xi1[accept]
# Set effect for 1st innings to be zero for idenitifiability
xi[1] <- 0
list(xi = xi, sigi = sigi)
} 

# Away effect parameter
awayUpdate <- function(data, nera, sigk, eta, lambda, theta, delta, agefun, xi, kappa, omega, phi, rwsd){
kappa1 <- rnorm(2, kappa, rwsd)
denrat <- (kappa^2 - kappa1^2)/(2*sigk) - (exp(-kappa1) - exp(-kappa))*with(data, tapply(eta[ID]*lambda*exp(-(theta[ID] + delta[YearScale] + agefun + xi[Inns] + omega[Hand] + phi[(Opposition - 1)*nera + Era])), HA, sum))
denrat <- denrat - (kappa1 - kappa)*with(data, tapply(eta[ID], HA, sum))
laccept <- pmin(0, denrat)
accept <- (log(runif(2)) < laccept)
kappa[accept] <- kappa1[accept]
# Set home effect to be zero for identifiability
kappa[1] <- 0
list(kappa = kappa, sigk = sigk)
}

# Handedness effect parameter
handUpdate <- function(data, nera, sigo, eta, lambda, theta, delta, agefun, xi, kappa, omega, phi, rwsd){
# omega1 <- rnorm(2, omega, rwsd)
# denrat <- (omega^2 - omega1^2)/(2*sigo) - (exp(-omega1) - exp(-omega))*with(data, tapply(eta[ID]*lambda*exp(-(theta[ID] + delta[YearScale] + agefun + xi[Inns] + kappa[HA] + phi[(Opposition - 1)*nera + Era])), Hand, sum))
# denrat <- denrat - (omega1 - omega)*with(data, tapply(eta[ID], Hand, sum))
# laccept <- pmin(0, denrat)
# accept <- (log(runif(2)) < laccept)
# omega[accept] <- omega1[accept]
# # Set right-hand to be zero for identifiability
omega[1:2] <- 0
list(omega = omega, sigo = sigo)
}

# Opposition effects parameters - now each team in each era
oppUpdate <- function(data, nopp, nera, oppdataInd, sigopp, eta, lambda, theta, delta, agefun, xi, kappa, omega, phi, rwsd){
phi1 <- rnorm(nopp*nera, phi, rwsd)
print(length(phi1))
denrat <- (phi^2 - phi1^2)/(2*sigopp)
denrat[oppdataInd] <- denrat[oppdataInd] - (exp(-phi1[oppdataInd]) - exp(-phi[oppdataInd]))*with(data, aggregate(eta[ID]*lambda*exp(-(theta[ID] + delta[YearScale] + agefun + xi[Inns] + kappa[HA] + omega[Hand])), list(Era, Opposition), sum))[, 3] - (phi1[oppdataInd] - phi[oppdataInd])*with(data, aggregate(eta[ID], list(Era, Opposition), sum))[, 3]
laccept <- pmin(0, denrat)
accept <- (log(runif(nopp)) < laccept)
phi[accept] <- phi1[accept]
# Impose corner constraint on Australia in most recent era
phi[nera] <- 0 
# phi[1:(nopp*nera)] <- 0 # To set all as zero
list(phi = phi, sigopp = sigopp)
} 

# Zero-inflation probability update
pizipUpdate <- function(data, np, zeroInd, D, censInd, pinns, pducks, pizipa, pizipb, lambda, pizip, rwsd){
pizip1logit <- rnorm(np, log(pizip/(1 - pizip)), rwsd)
pizip1 <- exp(pizip1logit)/(1 + exp(pizip1logit))
Ni <- pinns
Di <- pducks
denrat <- pizipa*log(pizip1/pizip) + (pizipb - 1 + Ni - Di)*log((1 - pizip1)/(1 - pizip)) + with(data, tapply(zeroInd*log(pizip1[ID] + (1 - pizip1[ID])*exp(-lambda*(1 - censInd))), ID, sum)) - with(data, tapply(zeroInd*log(pizip[ID] + (1 - pizip[ID])*exp(-lambda*(1 - censInd))), ID, sum)) + log(1 - pizip1) - log(1 - pizip)
laccept <- pmin(0, denrat)
accept <- (log(runif(np)) < laccept)
pizip[accept] <- pizip1[accept]
# print(summary(pizip))
# pizip <- rep(0.05, np)
list(pizip = pizip, pizipa = pizipa, pizipb = pizipb)
}

# Get log-likelihood as output
lhood <- function(data, nera, theta, delta, alphavec, xi, kappa, omega, phi, lambda, censInd, zeroInd){
poismu2 <- with(data, exp(theta[ID] + delta[YearScale] + alphavec[ID,] + xi[Inns] + kappa[HA] + omega[Hand] + phi[(Opposition - 1)*nera + Era]))
#llfull <- with(data, log(lambda*Inns) - lambda*Inns*Runs - sum(log(sequence(data$Runs)))
# Final term is a constant so ignore (can append afterwards if required)
# The likelihood needs censoring and ZIP components
ll2 <- with(data, sum(Runs*log(poismu2) - poismu2))
ll3 <- with(data, sum(Runs*log(lambda) - lambda))
#print(c(ll2, ll3))
list(loglik = ll2, loglik2 = ll3)
}
