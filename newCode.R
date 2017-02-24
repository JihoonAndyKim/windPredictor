


#load.R
rm(list=ls())
library(ncdf4)
library(abind)

options(scipen=999)

#Format the observations
formatDate <- function(obs) {
  #Time conversion
  obs$Date <- as.POSIXct(strptime(paste(obs$Date, obs$Hour, sep = " "), tz ="America/Denver", "%m/%d/%Y %H"))
  obs$Date <- format(obs$Date, tz="UTC", usetz = FALSE)
  obs$Date <- gsub("-| |:|[^0-9]", "", obs$Date)
  obs$Date <- substr(obs$Date, 1, nchar(obs$Date) - 4)
  
  obs <- obs[as.numeric(obs$Date) %% 100 == 0,]
  obs$Wind.Speed <- as.numeric(as.character(obs$Wind.Speed))
  return(obs)
}

u <- nc_open("u80mOneThree.nc")
v <- nc_open("v80mOneThree.nc")

u_time <- ncvar_get(u,"intTime")
u_valid <- ncvar_get(u, "intValidTime")
u_data <- ncvar_get(u,"U-component_of_wind_height_above_ground_80m")
v_time <- ncvar_get(v,"intTime")
v_valid <- ncvar_get(v, "intValidTime")
v_data <- ncvar_get(v,"V-component_of_wind_height_above_ground_80m")

#Get the absolute wind speed
absoluteWind <- sqrt(u_data^2 + v_data^2)
directionWind <- atan2(u_data, v_data)*180/pi
absoluteWind <- abind(aperm(absoluteWind, c(1,3,2)), v_valid)
directionWind <- abind(aperm(directionWind, c(1,3,2)), v_valid)

#Get the observed windspeeds
obs <- formatDate(read.csv("WindData.csv"))

#find just the 24 hour, 48 hour, and 72 hour forecasts
testWind <- absoluteWind[c(1,9,17),,]
testDirWind <- directionWind[c(1,9,17),,]

#binding an empty array to put the obs in
emptyBind = array(rep(0, 3 * 7232 * 2, 1), dim=c(3, 7232, 2))
testWind <- abind(testWind, emptyBind)
testDirWind <- abind(testDirWind, emptyBind)

#Replace all the time stamps with the appropriate observations and push the timestamps back one cell
for(i in seq(length(testWind[1,,12])))
{
  for(j in 1:3) {
    if(length(which(obs$Date == testWind[j,i,12])) != 0) {
      temp <- testWind[j, i, 12]
      testWind[j,i,12] <- obs[which(obs$Date == testWind[j,i,12]), 3]
      testWind[j,i,13] <- temp
    }
    else {
      temp <- testWind[j, i, 12]
      testWind[j,i,12] <- NA
      testWind[j,i,13] <- temp
    }
  }
}

for(i in seq(length(testDirWind[1,,12])))
{
  for(j in 1:3) {
    if(length(which(obs$Date == testDirWind[j,i,12])) != 0) {
      temp <- testDirWind[j,i,12]
      testDirWind[j,i,12] <- obs[which(obs$Date == testDirWind[j,i,12]), 3]
      testDirWind[j,i,13] <- temp
    }
    else {
      temp <- testDirWind[j,i,12]
      testDirWind[j,i,12] <- NA
      testDirWind[j,i,13] <- temp
    }
  }
}

#Put in seasons in the data
for(fore in 1:3) {
  for(days in 1:7232) {
    #extract month
    
    mon <- as.numeric(testWind[fore, days, 13]) %/% 10000 %% 100
    if(mon >= 12 || mon <= 2) {
      testWind[fore, days, 14] <- "WINTER"
    }
    else if(mon >= 3 && mon <= 5){
      testWind[fore, days, 14] <- "SPRING"
    }
    else if(mon >= 6 && mon <= 8) {
      testWind[fore, days, 14] <- "SUMMER"
    }
    else{
      testWind[fore, days, 14] <- "FALL"
    }
  }
}

for(fore in 1:3) {
  for(days in 1:7232) {
    #extract month
    mon <- as.numeric(testDirWind[fore, days, 13]) %/% 10000 %% 100
    if(mon >= 12 || mon <= 2) {
      testDirWind[fore, days, 14] <- "WINTER"
    }
    else if(mon >= 3 && mon <= 5){
      testDirWind[fore, days, 14] <- "SPRING"
    }
    else if(mon >= 6 && mon <= 8) {
      testDirWind[fore, days, 14] <- "SUMMER"
    }
    else{
      testDirWind[fore, days, 14] <- "FALL"
    }
  }
}

#Format it so that the gefs.R can read it
testWind <- aperm(testWind, c(2,3,1))
testDirWind <- aperm(testDirWind, c(2,3,1))

#Make all the values in our dataset numeric
for(i in 1:7232) {
  for(j in 1:13) {
    for(k in 1:3) {
      testWind[i,j,k] <- as.numeric(testWind[i,j,k])
      testDirWind[i,j,k] <- as.numeric(testDirWind[i,j,k])
    }
  }
}

save(testWind, file = "GEFS_NREL_Wind_Data.RData")
save(testDirWind, file = "GEFS_NREL_Wind_Dir_Data.RData")
nc_close(v)
rm(u)
rm(v)

#is(mydata)
#dim(mydata)











#gefs.R
#source('load.R')

removeNA <- function (testWind) {
  
  flag <- vector()
  for(i in seq(nrow(testWind[,,1]))) {
    if(all(is.na(testWind[i,12,]))) {
      flag <- c(flag, i)
    }
  }
  return(testWind[-flag,,])
  
}

windData = testWind

crps.gamma <- function(beta, obs, ensMean, ensVar){
  use <- !is.na(obs) & !is.na(ensMean)
  
  mu <- beta[1] + beta[2]*ensMean[use]
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


crps.truncnorm <- function(beta, obs, ensMean, ensVar){
  
  use <- !is.na(obs) & !is.na(ensMean)
  
  mu <- beta[1] + beta[2]*ensMean[use]
  sigma <- beta[3] + beta[4]*ensVar[use]
  
  y.tld <- (obs[use]-mu)/sqrt(sigma)
  c.tld <- mu/sqrt(sigma)
  
  crps <- sqrt(sigma)/pnorm(c.tld)^2 * (y.tld*pnorm(c.tld)*(2*pnorm(y.tld)+pnorm(c.tld)-2) +
                                          2*dnorm(y.tld)*pnorm(c.tld) - pnorm(sqrt(2)*c.tld)/sqrt(pi)  )
  
  return(mean(crps))
}

init.dates <- as.numeric(obs$Date)
months <- unique(init.dates %/% 10000)
nm <- length(months)

#Parameters as follows: month, location, lead time, parameters
# optVar <- array(dim=c(nm-36,5,3,4))
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
#       par0 <- c(0.1, 1, 0.1, 1)
#       mean2 <- function(x) {
#         mean(as.numeric(x))
#       }
#       var2 <- function(x) {
#         var(as.numeric(x))
#       }
#       
#       est <- optim(par0, crps.gamma,
#                    obs       = windData[ind.train,12,lt],
#                    ensMean   = apply(windData[ind.train,1:11,lt],1,mean2),
#                    ensVar    = apply(windData[ind.train,1:11,lt],1,var2),
#                    method    = "L-BFGS-B",
#                    lower     = c(2e-3, 0, 2e-3, 0),
#                    upper     = c(  15, 5,   20, 5),
#                    control = list(factr=1e-5/.Machine$double.eps))
#       
#       par0 <- est$par
#       optVar[im,st,lt,] <- est$par
#     }
#   }
# }

lowerBound = 0.4
upperBound = 0.7

startWindow = 5900
endWindow = 6000
count <- 0

for(testDay in startWindow:endWindow) {
  
  findTrainMonth = match(init.dates[testDay] %/% 10000, months)
  
  for(leadTime in 1:3) {
    
    testPar <- optVar[findTrainMonth - 36,1,leadTime,]
    testForecast <- windData[testDay,1:11,leadTime]
    mu <-testPar[1] * testPar[2]*mean(as.numeric(testForecast))
    sigma <- testPar[3] + testPar[4]*var(as.numeric(testForecast))
    
    lambda = mu/sigma
    alpha = mu^2/sigma
    
    variance = alpha/lambda^2
    print(variance)
    
    observed = windData[testDay + 165, 12, leadTime]
    
    percent <- qgamma(c(lowerBound, upperBound), shape=alpha, rate=lambda)
    
    # curve(dgamma(x, shape=alpha, rate=lambda), from=0, to=15, n=500, type="l", xlab = "Wind Speed", ylab = "")
    # abline(v = observed)
    # abline(v = percent[1], col = 'blue')
    # abline(v = percent[2], col = 'blue')
    
    if(!all(is.na(percent))) {
      if(observed > percent[1] & observed < percent[2]) {
        count = count + 1
      }
      else {
        # curve(dgamma(x, shape=alpha, rate=lambda), from=0, to=15, n=500, type="l", xlab = "Wind Speed", ylab = "")
        # abline(v = observed)
        # abline(v = percent[1], col = 'blue')
        # abline(v = percent[2], col = 'blue')
      }
    }
  }
}
print(count/300)
# plot(diffs, scores, type="o")
# abline(h = upperBound - lowerBound, col = 'blue')