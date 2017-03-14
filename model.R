#source('load.R')
library(circular)

removeNA <- function (testWind) {
  
  flag <- vector()
  for(i in seq(nrow(testWind[,,1]))) {
    if(all(is.na(testWind[i,12,]))) {
      flag <- c(flag, i)
    }
  }
  return(testWind[-flag,,])
  
}
mae <- function(error)
{
  mean(abs(error))
}
rmse <- function(error)
{
  sqrt(mean(error^2))
}


windData = testWind
windDir  = testDirWind

crps.gamma <- function(beta, obs, dirMean, ensMean, ensVar){
  use <- !is.na(obs) & !is.na(ensMean) & !is.na(dirMean)
  
  dirMean <- dirMean * pi / 180
  
  mu <- beta[1] + beta[2]*ensMean[use] + beta[5] * dvonmises(circular(dirMean[use]), mu = circular(beta[6]), kappa = beta[7])
  sigma <- beta[3] + beta[4]*ensVar[use]
  
  shape <- mu^2/sigma
  scale <- sigma/mu
  y <- as.numeric(obs[use])
  #print("------------------------")
  
  crps1 <- y*(2*pgamma(y,shape=shape,scale=scale)-1) - shape*scale*(2*pgamma(y,shape=shape+1,scale=scale)-1)
  crps2 <- shape*scale*beta(.5,shape+.5)/pi
  crps <- crps1 - crps2
  
  return(mean(crps))
}


# crps.truncnorm <- function(beta, obs, ensMean, ensVar){
#   
#   use <- !is.na(obs) & !is.na(ensMean)
#   
#   mu <- beta[1] + beta[2]*ensMean[use]
#   sigma <- beta[3] + beta[4]*ensVar[use]
#   
#   y.tld <- (obs[use]-mu)/sqrt(sigma)
#   c.tld <- mu/sqrt(sigma)
#   
#   crps <- sqrt(sigma)/pnorm(c.tld)^2 * (y.tld*pnorm(c.tld)*(2*pnorm(y.tld)+pnorm(c.tld)-2) +
#                                           2*dnorm(y.tld)*pnorm(c.tld) - pnorm(sqrt(2)*c.tld)/sqrt(pi)  )
#   
#   return(mean(crps))
# }








# init.dates <- as.numeric(obs$Date)
# months <- unique(init.dates %/% 10000)
# nm <- length(months)
# 
# #Parameters as follows: month, location, lead time, parameters
# optVar <- array(dim=c(nm-36,5,3,7))
# 
# for (st in 1:1)  {
#   cat(paste("Fitting station", st, "\n"))
#   for (lt in 1:3)  {
#     cat(paste("Fitting lead time", lt, "\n"))
#     for (im in 1:(nm-36))  {
#       print(im)
#       yyyy <- months[im+36] %/% 100
#       mm <- months[im+36] %% 100
#       yeardiff <- yyyy - months %/% 100
#       monthdiff <- pmin( abs(mm-months%%100), abs(mm-12-months%%100), abs(mm+12-months%%100) )
#       months.train <- months[yeardiff>0 & yeardiff<4 & monthdiff<=1]
#       ind.train <- (init.dates %/% 10000) %in% months.train
#       par0 <- c(0.1, 1, 0.1, 1, 0.1, 0, 1)
#       mean2 <- function(x) {
#         mean(as.numeric(x))
#       }
#       var2 <- function(x) {
#         var(as.numeric(x))
#       }
# 
#       est <- optim(par0, crps.gamma,
#                    obs       = windData[ind.train,12,lt],
#                    dirMean   = apply(windDir[ind.train,1:11,lt],1,mean2),
#                    ensMean   = apply(windData[ind.train,1:11,lt],1,mean2),
#                    ensVar    = apply(windData[ind.train,1:11,lt],1,var2),
#                    method    = "L-BFGS-B",
#                    lower     = c(2e-3, 0, 2e-3, 0, -pi, 0),
#                    upper     = c(  15, 5,   20, 5, 20, pi, 100),
#                    control = list(factr=1e-5/.Machine$double.eps))
# 
#       par0 <- est$par
#       optVar[im,st,lt,] <- est$par
#     }
#   }
# }
# 








lowerBound = 0.05
upperBound = 0.95

startWindow = 5900
endWindow = 6000
count <- 0
intervals <- c()
storeCRPSvonMises <- c()
Error <- 0
for(testDay in startWindow:endWindow) {
  
  findTrainMonth = match(init.dates[testDay] %/% 10000, months)
  
  #TODO: Change lead time to 3 or so
  for(leadTime in 1:1) {
    
    testPar <- optVar[findTrainMonth - 36,1,leadTime,]
    testForecast <- windData[testDay,1:11,leadTime]
    testWindDir <- mean(as.numeric(windDir[testDay, 1:11, leadTime])) * pi / 180
    
    d = dvonmises(circular(testWindDir), mu = circular(testPar[6]), kappa = testPar[7])
    mu <-testPar[1] * testPar[2]*mean(as.numeric(testForecast)) + testPar[5] * d
    sigma <- testPar[3] + testPar[4]*var(as.numeric(testForecast))
    
    lambda = mu/sigma
    alpha = mu^2/sigma
    
    shape = alpha
    scale = 1/lambda
    
    observed = as.numeric(windData[testDay + 165, 12, leadTime])
    meanOfGamma <- shape*scale
    
    Error <- Error + rmse(meanOfGamma - observed)
      
    crps1 <- observed*(2*pgamma(observed,shape=shape,scale=scale)-1) - shape*scale*(2*pgamma(observed,shape=shape+1,scale=scale)-1)
    crps2 <- shape*scale*beta(.5,shape+.5)/pi
    crps <- crps1 - crps2
    
    storeCRPSvonMises <- c(storeCRPSvonMises, crps)
    
    percent <- qgamma(c(lowerBound, upperBound), shape=alpha, rate=lambda)
    intervalWidth <- percent[2] - percent[1]
    intervals <- c(intervals, intervalWidth)
    
    if(!all(is.na(percent))) {
      if(observed > percent[1] & observed < percent[2]) {
        count = count + 1
      }
    }
  }
}
print(count/100)
print(Error)
print(mean(intervals))
#gamma interval widths and coverages
#MAE (median of gamma - observation) if not MSE (obs versus mean of gamma)

#Coverage from 0.05 to 0.95 theoretical is 0.87
#MAE = 257.8614


#Check to see if there is one bump
#Compute the bias in different directions
#Exploratory plots
#Interaction terms between von mises dist and forecast mean
#Refit the model without the intercept term (does it decrease kappa)
#Interaction term instead of the forecast mean (beta1*vonmises + beta2(vonmises * forecast mean))
#Cross-validation to other data
#use more than just four years
#Take a look to see if the means of the von mises are the same or no (force alphas to equal to each other)
#Use CRPS as well
#PIT histograms to check to see gamma and normal (to see which one is better)