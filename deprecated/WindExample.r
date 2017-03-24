load("./CUStuff/Research/Code/WindSpeed.Rdata")


crps.gamma <- function(param, obs, ensMean, ensVar){
  
  use <- !is.na(obs) & !is.na(ensMean)
  
  MEAN <- param[1] + param[2]*ensMean[use]
  VAR <- param[3] + param[4]*ensVar[use]
  
  shape <- MEAN^2/VAR
  scale <- VAR/MEAN
  y <- obs[use]
  
  crps1 <- y*(2*pgamma(y,shape=shape,scale=scale)-1) - shape*scale*(2*pgamma(y,shape=shape+1,scale=scale)-1)
  crps2 <- shape*scale*beta(.5,shape+.5)/pi
  crps <- crps1 - crps2

  return( mean(crps) )
}


crps.truncnorm <- function(param, obs, ensMean, ensVar){
  
  use <- !is.na(obs) & !is.na(ensMean)
  
  MEAN <- param[1] + param[2]*ensMean[use]
  VAR <- param[3] + param[4]*ensVar[use]
  
  y.tld <- (obs[use]-MEAN)/sqrt(VAR)
  c.tld <- MEAN/sqrt(VAR)
  
  crps <- sqrt(VAR)/pnorm(c.tld)^2 * (y.tld*pnorm(c.tld)*(2*pnorm(y.tld)+pnorm(c.tld)-2) +
             2*dnorm(y.tld)*pnorm(c.tld) - pnorm(sqrt(2)*c.tld)/sqrt(pi)  )

  return( mean(crps) )
}




# Model fitting

init.dates <- as.numeric(dimnames(wind.data[[1]])[[1]])
months <- unique(init.dates %/% 10000)
nm <- length(months)

par.EMOS <- array(dim=c(nm-36,5,3,4))

for (st in 1:5)  {
	cat(paste("Fitting station", st, "\n"))
	for (lt in 1:3)  {
		cat(paste("Fitting lead time", lt, "\n"))
		for (im in 1:(nm-36))  {
			yyyy <- months[im+36] %/% 100
			mm <- months[im+36] %% 100
			yeardiff <- yyyy - months %/% 100
			monthdiff <- pmin( abs(mm-months%%100), abs(mm-12-months%%100), abs(mm+12-months%%100) )
			months.train <- months[yeardiff>0 & yeardiff<4 & monthdiff<=1]
			ind.train <- (init.dates %/% 10000) %in% months.train

			par0 <- c(0.1, 1, 0.1, 1)
  
    		est <- optim(par0, crps.truncnorm,  # crps.gamma,
                 obs       = wind.data[[st]][ind.train,12,lt],
                 ensMean   = apply(wind.data[[st]][ind.train,1:11,lt],1,mean,na.rm=TRUE),   
                 ensVar    = apply(wind.data[[st]][ind.train,1:11,lt],1,var,na.rm=TRUE),
                 method    = "L-BFGS-B",
                 lower     = c(2e-3, 0, 2e-3, 0),
                 upper     = c(  15, 5,   20, 5),
                 control = list(factr=1e-5/.Machine$double.eps))

			par0 <- est$par
			par.EMOS[im,st,lt,] <- est$par
		}
	}
}

save(par.EMOS, file="EMOS-Parameters-TrNorm.Rdata")




# Calculate PIT values

verif.ind <- (init.dates%/%1000000) > 2002
l <- sum(verif.ind)

pit.EMOS <- array(dim=c(l,5,3))


for (st in 1:5)  {
	cat(paste("Fitting station", st, "\n"))
	for (lt in 1:3)  {
		cat(paste("Fitting lead time", lt, "\n"))
		for (im in 1:(nm-36))  {
			ind.current.month <- (init.dates %/% 10000) %in% months[im+36]

         obs     <- wind.data[[st]][ind.current.month,12,lt]
         ensMean <- apply(wind.data[[st]][ind.current.month,1:11,lt],1,mean,na.rm=TRUE)   
         ensVar  <- apply(wind.data[[st]][ind.current.month,1:11,lt],1,var,na.rm=TRUE)

			pred.mean <- par.EMOS[im,st,lt,1] + par.EMOS[im,st,lt,2]*ensMean
			pred.var  <- par.EMOS[im,st,lt,3] + par.EMOS[im,st,lt,4]*ensVar

#			shape <- pred.mean^2/pred.var
#			scale <- pred.var/pred.mean

			ind.match <- ind.current.month[verif.ind]
#			pit.EMOS[ind.match,st,lt] <- pgamma(obs,shape=shape,scale=scale)
			pit.EMOS[ind.match,st,lt] <- (pnorm(obs,pred.mean,sqrt(pred.var))
						-pnorm(0,pred.mean,sqrt(pred.var))) / pnorm(pred.mean/sqrt(pred.var))
		}
	}
}



K <- 13

par(mfrow=c(3,5))

for (lt in 1:3)  {
	for (st in 1:5)  {
		hist(pit.EMOS[,st,lt], breaks=(0:(K+1))/(K+1), xaxt="n",
			main=paste("Wind park ",st,", ", 24*lt, "h", sep=""), col="lightblue",
			cex.main=1.5, cex.lab=1.2, xlab="probability integral transform", ylab="frequency")
			axis(1, at=seq(0, 1, 0.2), labels=(0:5)/5)
    		abline(h=sum(!is.na(pit.EMOS[,st,lt]))/(K+1), lty=2)
	}
}



##
#   Now: multivariate verification 
#    (ensemble vs. ECC-Q vs. Schaake-Q vs. ECC-R110 vs. Schaake-R110)
 

gini.md <- function(x,na.rm=FALSE)  {
   if(na.rm & any(is.na(x)))  x <- x[!is.na(x)] 
   n <-length(x)
   return(4*sum((1:n)*sort(x,na.last=TRUE))/(n^2)-2*mean(x)*(n+1)/n)             # biased version!
}

crps.score <- function(obs,fcst)  {
	return( mean(abs(fcst-obs)) - gini.md(fcst)/2 )
}


variogram.score <- function(obs,fcst,wgt,p=2)  {
	obs.diff <- abs(as.vector(dist(obs)))^p
	fcst.diff <- apply(abs(apply(fcst,2,function(x) as.vector(dist(x))))^p, 1, mean)
	return( sum(wgt*(obs.diff-fcst.diff)^2) )
}

energy.score <- function(obs,fcst)  {
	T1 <- mean(sqrt(apply(sweep(fcst,1,obs,"-")^2,2,mean)))
	T2 <- mean(sqrt(apply(apply(fcst,1,diff)^2,1,mean)))
	return(T1-0.5*T2)
}


qtruncnorm <- function(p, mu, sigma2)  {
	sigma <- sqrt(sigma2)
	return( mu + sigma*qnorm(1+pnorm(mu/sigma)*(p-1)) )
}


lon.stat <- c(-103.662,-104.042,-103.388,-102.636,-102.307)
lat.stat <- c(  39.373,  40.853,  40.994,  37.700,  39.010)


gcd.slc <- function(coords1, coords2=coords1) {
   R <- 6371 # Earth mean radius [km]
   deg2rad <- pi/180
   lon1 <- coords1[,1]*deg2rad
   lat1 <- coords1[,2]*deg2rad
   lon2 <- coords2[,1]*deg2rad
   lat2 <- coords1[,2]*deg2rad
   d <- acos( outer(sin(lat1),sin(lat2),"*") +
		 outer(cos(lat1),cos(lat2),"*") * cos(outer(lon2,lon1,"-")) ) * R
   d[outer(lon1,lon2,"==") & outer(lat1,lat2,"==")] <- 0
   return(d) # Distance in km
}


method.names <- c("ensemble","ECC-Q","Schaake-Q","Random-Q","Ordered-Q")
ES <- VS.05 <- VS.1 <- array(dim=c(l,5,3), dimnames=list(init.dates[verif.ind],method.names,(1:3)*24))
CRPS.min <- CRPS.max <- array(dim=c(l,5,3), dimnames=list(init.dates[verif.ind],method.names,(1:3)*24))

dst <- gcd.slc(cbind(lon.stat,lat.stat))
wgt <- 1/dst[lower.tri(dst)]

# sample.days <- matrix(NA,10,11)


for (im in 1:(nm-36))
{
	cat(paste("Calculating scores for ", months[im+36]%%100, "/", months[im+36]%/%100, "\n", sep=""))

#	yyyy <- months[im+36] %/% 100
#	mm <- months[im+36] %% 100
#	yeardiff <- yyyy - months %/% 100
#	monthdiff <- pmin( abs(mm-months%%100), abs(mm-12-months%%100), abs(mm+12-months%%100) )
#	months.train <- months[yeardiff>0 & yeardiff<4 & monthdiff<=1]
#	ind.train <- (init.dates %/% 10000) %in% months.train

	ind.current.month <- (init.dates %/% 10000) %in% months[im+36]
	pred.mean <- pred.var <- verif.obs <- matrix(NA,sum(ind.current.month),5)
	raw.ens <- array(dim=c(sum(ind.current.month),11,5))

	ind.match <- which(ind.current.month[verif.ind])
	schaake.ens <- ecc.ens <- random.ens <- ordered.ens <- matrix(NA,11,5)
#	schaake.r110.ens <- ecc.r110.ens <- matrix(NA,110,5)

	for (lt in 1:3)  {
		for (st in 1:5)  {
			ensMean <- apply(wind.data[[st]][ind.current.month,1:11,lt],1,mean,na.rm=TRUE)
			ensVar  <- apply(wind.data[[st]][ind.current.month,1:11,lt],1,var,na.rm=TRUE)
			pred.mean[,st]  <- par.EMOS[im,st,lt,1] + par.EMOS[im,st,lt,2]*ensMean
			pred.var[,st]   <- par.EMOS[im,st,lt,3] + par.EMOS[im,st,lt,4]*ensVar
			verif.obs[,st]  <- wind.data[[st]][ind.current.month,12,lt]
			raw.ens[,,st]   <- wind.data[[st]][ind.current.month,1:11,lt]
		}

		for (dd in seq_along(ind.match))  {
		  # ensemble
			ES[ind.match[dd],1,lt]    <- energy.score(verif.obs[dd,],t(raw.ens[dd,,]))
			VS.05[ind.match[dd],1,lt] <- variogram.score(verif.obs[dd,],t(raw.ens[dd,,]),wgt,0.5)
			VS.1[ind.match[dd],1,lt]  <- variogram.score(verif.obs[dd,],t(raw.ens[dd,,]),wgt,1)

			CRPS.min[ind.match[dd],1,lt] <- crps.score(min(verif.obs[dd,]),apply(raw.ens[dd,,],1,min))
			CRPS.max[ind.match[dd],1,lt] <- crps.score(max(verif.obs[dd,]),apply(raw.ens[dd,,],1,max))


		  # ECC-Q, Schaake-Q, Random-Q, Ordered-Q
			recent.ind <- tail(which(init.dates<init.dates[ind.current.month][dd]),11)
   		for (st in 1:5)  {
				recent.obs <- wind.data[[st]][recent.ind,12,lt]
				EMOS.quantiles <- qtruncnorm(seq(1,11,1)/12,pred.mean[dd,st],pred.var[dd,st])
		      ecc.ens[order(raw.ens[dd,,st]),st] <- EMOS.quantiles
				schaake.ens[order(recent.obs),st] <- EMOS.quantiles
				random.ens[sample(1:11,11),st] <- EMOS.quantiles
				ordered.ens[,st] <- EMOS.quantiles
  			}
			ES[ind.match[dd],2,lt]    <- energy.score(verif.obs[dd,],t(ecc.ens))
			VS.05[ind.match[dd],2,lt] <- variogram.score(verif.obs[dd,],t(ecc.ens),wgt,0.5)
			VS.1[ind.match[dd],2,lt]  <- variogram.score(verif.obs[dd,],t(ecc.ens),wgt,1)

			CRPS.min[ind.match[dd],2,lt] <- crps.score(min(verif.obs[dd,]),apply(ecc.ens,1,min))
			CRPS.max[ind.match[dd],2,lt] <- crps.score(max(verif.obs[dd,]),apply(ecc.ens,1,max))


			ES[ind.match[dd],3,lt]    <- energy.score(verif.obs[dd,],t(schaake.ens))
			VS.05[ind.match[dd],3,lt] <- variogram.score(verif.obs[dd,],t(schaake.ens),wgt,0.5)
			VS.1[ind.match[dd],3,lt]  <- variogram.score(verif.obs[dd,],t(schaake.ens),wgt,1)

			CRPS.min[ind.match[dd],3,lt] <- crps.score(min(verif.obs[dd,]),apply(schaake.ens,1,min))
			CRPS.max[ind.match[dd],3,lt] <- crps.score(max(verif.obs[dd,]),apply(schaake.ens,1,max))


			ES[ind.match[dd],4,lt]    <- energy.score(verif.obs[dd,],t(random.ens))
			VS.05[ind.match[dd],4,lt] <- variogram.score(verif.obs[dd,],t(random.ens),wgt,0.5)
			VS.1[ind.match[dd],4,lt]  <- variogram.score(verif.obs[dd,],t(random.ens),wgt,1)

			CRPS.min[ind.match[dd],4,lt] <- crps.score(min(verif.obs[dd,]),apply(random.ens,1,min))
			CRPS.max[ind.match[dd],4,lt] <- crps.score(max(verif.obs[dd,]),apply(random.ens,1,max))


			ES[ind.match[dd],5,lt]    <- energy.score(verif.obs[dd,],t(ordered.ens))
			VS.05[ind.match[dd],5,lt] <- variogram.score(verif.obs[dd,],t(ordered.ens),wgt,0.5)
			VS.1[ind.match[dd],5,lt]  <- variogram.score(verif.obs[dd,],t(ordered.ens),wgt,1)

			CRPS.min[ind.match[dd],5,lt] <- crps.score(min(verif.obs[dd,]),apply(ordered.ens,1,min))
			CRPS.max[ind.match[dd],5,lt] <- crps.score(max(verif.obs[dd,]),apply(ordered.ens,1,max))

#			for (rpt in 1:10)  {
#				sample.days[rpt,] <- sample(1:sum(ind.train),11)
#			}
#
#		  # ECC-R and Schaake-R
#   		for (st in 1:5)  {
#				hist.obs <- wind.data[[st]][ind.train,12,lt]
#				for (rpt in 1:10)  {
#					EMOS.sim <- qtruncnorm(runif(11),pred.mean[dd,st],pred.var[dd,st])
#		      	ecc.r110.ens[(rpt-1)*11+order(raw.ens[dd,,st]),st] <- EMOS.sim
#					schaake.r110.ens[(rpt-1)*11+order(hist.obs[sample.days[rpt,]]),st] <- EMOS.sim
#				}
# 			}
#			ES[ind.match[dd],4,lt]    <- energy.score(verif.obs[dd,],t(ecc.r110.ens))
#			VS.05[ind.match[dd],4,lt] <- variogram.score(verif.obs[dd,],t(ecc.r110.ens),wgt,0.5)
#			VS.1[ind.match[dd],4,lt]  <- variogram.score(verif.obs[dd,],t(ecc.r110.ens),wgt,1)
#
#			ES[ind.match[dd],5,lt]    <- energy.score(verif.obs[dd,],t(schaake.r110.ens))
#			VS.05[ind.match[dd],5,lt] <- variogram.score(verif.obs[dd,],t(schaake.r110.ens),wgt,0.5)
#			VS.1[ind.match[dd],5,lt]  <- variogram.score(verif.obs[dd,],t(schaake.r110.ens),wgt,1)

		}
	}
}


save(ES, VS.05, VS.1, CRPS.min, CRPS.max, file="EMOS-Scores.Rdata")




for (i in 1:3)  {
	dev.new()
	par(mfrow=c(1,5))
	boxplot(ES[,,i], col="lightblue")
	boxplot(VS.05[,,i], col="lightblue")
	boxplot(VS.1[,,i], col="lightblue")
	boxplot(CRPS.min[,,i], col="lightblue")
	boxplot(CRPS.max[,,i], col="lightblue")
}


ind <- tail(1:sum(verif.ind),365)

ES.ind <- apply(ES[ind,,],c(2,3),mean,na.rm=TRUE)
VS.05.ind <- apply(VS.05[ind,,],c(2,3),mean,na.rm=TRUE)
VS.1.ind <- apply(VS.1[ind,,],c(2,3),mean,na.rm=TRUE)


signif(t(ES.ind),3)

signif(t(VS.05.ind),3)

signif(t(VS.1.ind),3)


#    ensemble ECC-Q  Schaake-Q	Random-Q	 Ordered-Q
# 24   1.79    1.46    1.51       1.47      	2.29
# 48   1.78    1.57    1.62       1.59      	2.52
# 72   1.78    1.67    1.72       1.69      	2.66
#
#    ensemble  ECC-Q  Schaake-Q	Random-Q	 Ordered-Q
# 24  0.0265  0.0220   0.0225     0.0241    	0.0304
# 48  0.0259  0.0228   0.0236     0.0270    	0.0319
# 72  0.0253  0.0244   0.0251     0.0282    	0.0370
#
#    ensemble  ECC-Q  Schaake-Q	Random-Q	 Ordered-Q
# 24   0.252   0.214    0.218     0.226      0.267
# 48   0.251   0.226    0.230     0.256      0.287
# 72   0.248   0.242    0.244     0.266      0.322
#



# Skill scores (w.r.t. ensemble)

signif(1-ES.ind[-1,]/ES.ind[1,],3)

signif(1-VS.05.ind[-1,]/VS.05.ind[1,],3)

signif(1-VS.1.ind[-1,]/VS.1.ind[1,],3)

#               24     48     72
# ECC-Q      0.184  0.119  0.063
# Schaake-Q  0.152  0.092  0.037
# Random-Q   0.175  0.108  0.051
# Ordered-Q -0.284 -0.420 -0.493
#
# ECC-Q      0.171  0.119  0.036
# Schaake-Q  0.134  0.068  0.055
# Random-Q   0.047 -0.020 -0.087
# Ordered-Q -0.147 -0.231 -0.461
#
# ECC-Q      0.151  0.096  0.027
# Schaake-Q  0.131  0.073  0.031
# Random-Q   0.088 -0.017 -0.063
# Ordered-Q -0.062 -0.145 -0.299



# Significance tests

i <- 2

for (lt in 1:3)  {
   cat("\nlead time", lt, ifelse(lt==1,"day\n","days\n"))
	for (j in c(1,4,5))  {
		p.value <- round(wilcox.test(ES[ind,i,lt],ES[ind,j,lt],paired=TRUE)$p.value,4)
		cat(paste(dimnames(ES)[[2]][i], "vs.", dimnames(ES)[[2]][j], ":", p.value,"\n"))
		delta <- ES[ind,i,lt] - ES[ind,j,lt]
		cat(paste("lag 1 correlation: ", cor(rank(delta[-365]),rank(delta[-1])),"\n\n"))
	}
}

for (lt in 1:3)  {
   cat("\nlead time", lt, ifelse(lt==1,"day\n","days\n"))
	for (j in c(1,4,5))  {
		p.value <- round(wilcox.test(VS.05[ind,i,lt],VS.05[ind,j,lt],paired=TRUE)$p.value,4)
		cat(paste(dimnames(ES)[[2]][i], "vs.", dimnames(ES)[[2]][j], ":", p.value,"\n"))
		delta <- VS.05[ind,i,lt] - VS.05[ind,j,lt]
		cat(paste("lag 1 correlation: ", cor(rank(delta[-365]),rank(delta[-1])),"\n\n"))
	}
}

for (lt in 1:3)  {
   cat("\nlead time", lt, ifelse(lt==1,"day\n","days\n"))
	for (j in c(1,4,5))  {
		p.value <- round(wilcox.test(VS.1[ind,i,lt],VS.1[ind,j,lt],paired=TRUE)$p.value,4)
		cat(paste(dimnames(ES)[[2]][i], "vs.", dimnames(ES)[[2]][j], ":", p.value,"\n"))
		delta <- VS.1[ind,i,lt] - VS.1[ind,j,lt]
		cat(paste("lag 1 correlation: ", cor(rank(delta[-365]),rank(delta[-1])),"\n\n"))
	}
}



par(mfrow=c(3,1))
plot(ES[ind,2,1]-ES[ind,4,1], pch=19, cex=0.3)
abline(h=0,col=2)
plot(VS.05[ind,2,1]-VS.05[ind,4,1], pch=19, cex=0.3)
abline(h=0,col=2)
plot(VS.1[ind,2,1]-VS.1[ind,4,1], pch=19, cex=0.3)
abline(h=0,col=2)































# Transform to Gaussian and estimate empirical copula

ptruncnorm <- function(x, mu, sigma2)  {
	sigma <- sqrt(sigma2)
	return( (pnorm(x,mu,sigma)-pnorm(0,mu,sigma)) / pnorm(mu/sigma) )
}


GC.matrix <- array(dim=c(15,15,nm-36))

for (im in 1:(nm-36))  {
	yyyy <- months[im+36] %/% 100
	mm <- months[im+36] %% 100
	yeardiff <- yyyy - months %/% 100
	monthdiff <- pmin( abs(mm-months%%100), abs(mm-12-months%%100), abs(mm+12-months%%100) )
	months.train <- months[yeardiff>0 & yeardiff<4 & monthdiff<=1]
	ind.train <- (init.dates %/% 10000) %in% months.train

	cat(paste("Fitting copula for ", mm, "/", yyyy, "\n", sep=""))

	error.stdnorm <- matrix(NA,sum(ind.train),15)

	for (st in 1:5)  {
		for (lt in 1:3)  {
			obs     = wind.data[[st]][ind.train,12,lt]
			ensMean = apply(wind.data[[st]][ind.train,1:11,lt],1,mean,na.rm=TRUE)
			ensVar  = apply(wind.data[[st]][ind.train,1:11,lt],1,var,na.rm=TRUE)

			pred.mean <- par.EMOS[im,st,lt,1] + par.EMOS[im,st,lt,2]*ensMean
			pred.var  <- par.EMOS[im,st,lt,3] + par.EMOS[im,st,lt,4]*ensVar

			error.stdnorm[,(st-1)*3+lt] <- qnorm(ptruncnorm(obs, pred.mean, pred.var))
		}
	}

	GC.matrix[,,im] <- crossprod(error.stdnorm,error.stdnorm) / sum(ind.train)
	image(GC.matrix[,,im], col=tim.colors(50))
}

save(GC.matrix, file="CopulaMatrices-TrNorm.Rdata")






### Now: multivariate verification 


variogram.score <- function(obs,fcst,wgt,p=2)  {
	obs.diff <- abs(as.vector(dist(obs)))^p
	fcst.diff <- apply(abs(apply(fcst,2,function(x) as.vector(dist(x))))^p, 1, mean)
	return( sum(wgt*(obs.diff-fcst.diff)^2) )
}

energy.score <- function(obs,fcst)  {
	T1 <- mean(sqrt(apply(sweep(fcst,1,obs,"-")^2,2,mean)))
	T2 <- mean(sqrt(apply(apply(fcst,1,diff)^2,1,mean)))
	return(T1-0.5*T2)
}


qtruncnorm <- function(p, mu, sigma2)  {
	sigma <- sqrt(sigma2)
	return( mu + sigma*qnorm(1+pnorm(mu/sigma)*(p-1)) )
}


lon.stat <- c(-103.662,-104.042,-103.388,-102.636,-102.307)
lat.stat <- c(  39.373,  40.853,  40.994,  37.700,  39.010)


gcd.slc <- function(coords1, coords2=coords1) {
   R <- 6371 # Earth mean radius [km]
   deg2rad <- pi/180
   lon1 <- coords1[,1]*deg2rad
   lat1 <- coords1[,2]*deg2rad
   lon2 <- coords2[,1]*deg2rad
   lat2 <- coords1[,2]*deg2rad
   d <- acos( outer(sin(lat1),sin(lat2),"*") +
		 outer(cos(lat1),cos(lat2),"*") * cos(outer(lon2,lon1,"-")) ) * R
   d[outer(lon1,lon2,"==") & outer(lat1,lat2,"==")] <- 0
   return(d) # Distance in km
}


method.names <- c("ensemble","ECC-Q","GC-11","GC-1000")
ES <- VS.05 <- VS.1 <- array(dim=c(l,4,3), dimnames=list(init.dates[verif.ind],method.names,(1:3)*24))

dst <- gcd.slc(cbind(lon.stat,lat.stat))
wgt <- 1/dst[lower.tri(dst)]


for (im in 1:(nm-36))
{
	cat(paste("Calculating scores for ", months[im+36]%%100, "/", months[im+36]%/%100, "\n", sep=""))

	ind.current.month <- (init.dates %/% 10000) %in% months[im+36]
	pred.mean <- pred.var <- verif.obs <- matrix(NA,sum(ind.current.month),5)
	raw.ens <- array(dim=c(sum(ind.current.month),11,5))

	ind.match <- which(ind.current.month[verif.ind])
	ecc.ens <- matrix(NA,11,5)
	gc1000.ens <- matrix(NA,5,1000)

	for (lt in 1:3)  {
		for (st in 1:5)  {
			ensMean <- apply(wind.data[[st]][ind.current.month,1:11,lt],1,mean,na.rm=TRUE)
			ensVar  <- apply(wind.data[[st]][ind.current.month,1:11,lt],1,var,na.rm=TRUE)
			pred.mean[,st] <- par.EMOS[im,st,lt,1] + par.EMOS[im,st,lt,2]*ensMean
			pred.var[,st]  <- par.EMOS[im,st,lt,3] + par.EMOS[im,st,lt,4]*ensVar
			verif.obs[,st] <- wind.data[[st]][ind.current.month,12,lt]
			raw.ens[,,st]  <- wind.data[[st]][ind.current.month,1:11,lt]
		}

		for (dd in seq_along(ind.match))  {
		  # ensemble
			ES[ind.match[dd],1,lt]    <- energy.score(verif.obs[dd,],t(raw.ens[dd,,]))
			VS.05[ind.match[dd],1,lt] <- variogram.score(verif.obs[dd,],t(raw.ens[dd,,]),wgt,0.5)
			VS.1[ind.match[dd],1,lt]  <- variogram.score(verif.obs[dd,],t(raw.ens[dd,,]),wgt,1)

		  # ECC-Q
   		for (st in 1:5)  {
				EMOS.quantiles <- qtruncnorm(seq(0.5,10.5,1)/11,pred.mean[dd,st],pred.var[dd,st])
		      ecc.ens[order(raw.ens[dd,,st]),st] <- EMOS.quantiles
  			}
			ES[ind.match[dd],2,lt]    <- energy.score(verif.obs[dd,],t(ecc.ens))
			VS.05[ind.match[dd],2,lt] <- variogram.score(verif.obs[dd,],t(ecc.ens),wgt,0.5)
			VS.1[ind.match[dd],2,lt]  <- variogram.score(verif.obs[dd,],t(ecc.ens),wgt,1)

		  # GC-11 and GC-1000
  			R <- chol(cov2cor(GC.matrix[lt+(0:4)*3,lt+(0:4)*3,im]))
			Z <- crossprod(R,matrix(rnorm(5*1000),5,1000))
   		for (st in 1:5)  {
				gc1000.ens[st,] <- qtruncnorm(pnorm(Z[st,]),pred.mean[dd,st],pred.var[dd,st])
  		 	}
			ES[ind.match[dd],3,lt]    <- energy.score(verif.obs[dd,],gc1000.ens[,1:11])
			VS.05[ind.match[dd],3,lt] <- variogram.score(verif.obs[dd,],gc1000.ens[,1:11],wgt,0.5)
			VS.1[ind.match[dd],3,lt]  <- variogram.score(verif.obs[dd,],gc1000.ens[,1:11],wgt,1)

			ES[ind.match[dd],4,lt]    <- energy.score(verif.obs[dd,],gc1000.ens)
			VS.05[ind.match[dd],4,lt] <- variogram.score(verif.obs[dd,],gc1000.ens,wgt,0.5)
			VS.1[ind.match[dd],4,lt]  <- variogram.score(verif.obs[dd,],gc1000.ens,wgt,1)
		}
	}
}


save(ES, VS.05, VS.1, file="EMOS-Scores.Rdata")




for (i in 1:3)  {
	dev.new()
	par(mfrow=c(1,3))
	boxplot(ES[,,i], col="lightblue")
	boxplot(VS.05[,,i], col="lightblue")
	boxplot(VS.1[,,i], col="lightblue")
}


signif(t(apply(ES,c(2,3),mean,na.rm=TRUE)),3)

signif(t(apply(VS.05,c(2,3),mean,na.rm=TRUE)),3)

signif(t(apply(VS.1,c(2,3),mean,na.rm=TRUE)),3)

#    ensemble ECC-Q GC-11 GC-1000
# 24     1.89  1.37  1.48    1.48
# 48     1.92  1.51  1.62    1.62
# 72     1.89  1.61  1.70    1.71

#    ensemble  ECC-Q  GC-11 GC-1000
# 24   0.0270 0.0212 0.0218  0.0202
# 48   0.0268 0.0226 0.0233  0.0213
# 72   0.0256 0.0231 0.0237  0.0217

#    ensemble ECC-Q GC-11 GC-1000
# 24    0.245 0.198 0.206   0.191
# 48    0.250 0.213 0.224   0.204
# 72    0.241 0.218 0.228   0.208



# Significance test

lt <- 3

for (i in 1:3)  {
	for (j in (i+1):4)  {
		p.value <- round(wilcox.test(ES[,i,lt],ES[,j,lt],paired=TRUE)$p.value,4)
		cat(paste(dimnames(ES)[[2]][i], "vs.", dimnames(ES)[[2]][j], ":", p.value,"\n"))
	}
}

for (i in 1:3)  {
	for (j in (i+1):4)  {
		p.value <- round(wilcox.test(VS.05[,i,lt],VS.05[,j,lt],paired=TRUE)$p.value,4)
		cat(paste(dimnames(ES)[[2]][i], "vs.", dimnames(ES)[[2]][j], ":", p.value,"\n"))
	}
}

for (i in 1:3)  {
	for (j in (i+1):4)  {
		p.value <- round(wilcox.test(VS.1[,i,lt],VS.1[,j,lt],paired=TRUE)$p.value,4)
		cat(paste(dimnames(ES)[[2]][i], "vs.", dimnames(ES)[[2]][j], ":", p.value,"\n"))
	}
}




# Modified band depth histograms and average rank histograms


avg.rank <- function(x)  {
	 missing <- apply(is.na(x),1,any)
    x.ranks <- apply(x[!missing,],1,rank)
    x.preranks <- apply(x.ranks,1,mean)
    x.rank <- rank(x.preranks,ties="random")
    return(x.rank[1])
  }


method.names <- c("ensemble","ECC-Q","GC-11","GC-1000")
avg.rank.results <- array(dim=c(l,4,3), dimnames=list(init.dates[verif.ind],method.names,(1:3)*24))


for (im in 1:(nm-36))
{
	cat(paste("Calculating scores for ", months[im+36]%%100, "/", months[im+36]%/%100, "\n", sep=""))

	ind.current.month <- (init.dates %/% 10000) %in% months[im+36]
	pred.mean <- pred.var <- verif.obs <- matrix(NA,sum(ind.current.month),5)
	raw.ens <- array(dim=c(sum(ind.current.month),11,5))

	ind.match <- which(ind.current.month[verif.ind])
	ecc.ens <- matrix(NA,11,5)
	gc1000.ens <- matrix(NA,5,1000)

	for (lt in 1:3)  {
		for (st in 1:5)  {
			ensMean <- apply(wind.data[[st]][ind.current.month,1:11,lt],1,mean,na.rm=TRUE)
			ensVar  <- apply(wind.data[[st]][ind.current.month,1:11,lt],1,var,na.rm=TRUE)
			pred.mean[,st] <- par.EMOS[im,st,lt,1] + par.EMOS[im,st,lt,2]*ensMean
			pred.var[,st]  <- par.EMOS[im,st,lt,3] + par.EMOS[im,st,lt,4]*ensVar
			verif.obs[,st] <- wind.data[[st]][ind.current.month,12,lt]
			raw.ens[,,st]  <- wind.data[[st]][ind.current.month,1:11,lt]
		}

		for (dd in seq_along(ind.match))  {
		  # ensemble
			avg.rank.results[ind.match[dd],1,lt] <- avg.rank(cbind(verif.obs[dd,],t(raw.ens[dd,,])))

		  # ECC-Q
   		for (st in 1:5)  {
				EMOS.quantiles <- qtruncnorm(seq(0.5:10.5)/11,pred.mean[dd,st],pred.var[dd,st])
		      ecc.ens[order(raw.ens[dd,,st]),st] <- EMOS.quantiles
   		}
			avg.rank.results[ind.match[dd],2,lt] <- avg.rank(cbind(verif.obs[dd,],t(ecc.ens)))

		  # GC-11 and GC-1000
  			R <- chol(cov2cor(GC.matrix[lt+(0:4)*3,lt+(0:4)*3,im]))
			Z <- crossprod(R,matrix(rnorm(5*1000),5,1000))
   		for (st in 1:5)  {
				gc1000.ens[st,] <- qtruncnorm(pnorm(Z[st,]),pred.mean[dd,st],pred.var[dd,st])
  		 	}
			avg.rank.results[ind.match[dd],3,lt] <- avg.rank(cbind(verif.obs[dd,],gc1000.ens[,1:11]))
			avg.rank.results[ind.match[dd],4,lt] <- avg.rank(cbind(verif.obs[dd,],gc1000.ens))
		}
	}
}



for (lt in 1:3)  {
	dev.new()
	par(mfrow=c(1,4))
	for (j in 1:3)  {
      hist(avg.rank.results[,j,lt],xlim=c(0.5,12.5),breaks=seq(0.5,12.5,by=1),main=method.names[j],
         freq=FALSE,axes=FALSE,xlab="",ylab="",col="#4876FF40",border=FALSE)
      axis(1,at=c(1:12))
      axis(2)
      abline(h=1/12)
	}
   hist(avg.rank.results[,4,lt]/1002,xlim=c(0,1),breaks=seq(0,1,,12),main=method.names[j],
      freq=FALSE,axes=FALSE,xlab="",ylab="",col="#4876FF40",border=FALSE)
   axis(1,at=seq(1,23,2)/24, labels=(1:12))
   axis(2)
   abline(h=1)
}

