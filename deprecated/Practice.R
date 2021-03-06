load("WindSpeed.Rdata")

windData = wind.data[[1]]

crps.gamma <- function(beta, obs, ensMean, ensVar){
  
  use <- !is.na(obs) & !is.na(ensMean)
  
  mu <- beta[1] + beta[2]*ensMean[use]
  sigma <- beta[3] + beta[4]*ensVar[use]
  
  shape <- mu^2/sigma
  scale <- sigma/mu
  y <- obs[use]
  
  #print(y)
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

# init.dates <- as.numeric(dimnames(wind.data[[1]])[[1]])
# months <- unique(init.dates %/% 10000)
# nm <- length(months)
# 
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
# 
#       est <- optim(par0, crps.truncnorm,
#                    obs       = wind.data[[st]][ind.train,12,lt],
#                    ensMean   = apply(wind.data[[st]][ind.train,1:11,lt],1,mean,na.rm=TRUE),
#                    ensVar    = apply(wind.data[[st]][ind.train,1:11,lt],1,var,na.rm=TRUE),
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

lowerBound = 0.10
upperBound = 0.90

startWindow = 4000
scores <- vector()
diffs <- vector()

for(endWindow in  4010:4100) {
  diff = endWindow - startWindow + 1
  diffs = c(diffs, diff)
  count = 0
  
  for(testDay in startWindow:endWindow) {
    findTrainMonth = match(init.dates[testDay] %/% 10000, months)
    for(leadTime in 1:3) {
      testPar <- optVar[findTrainMonth - 36,1,leadTime,]
      testForecast <- windData[testDay,1:11,leadTime]
      mu <-testPar[1] * testPar[2]*mean(testForecast)
      sigma <- testPar[3] + testPar[4]*var(testForecast)
      
      lambda = mu/sigma
      alpha = mu^2/sigma
      
      observed = windData[testDay, 12, leadTime]
      
      percent <- qgamma(c(lowerBound, upperBound), shape=alpha, rate=lambda)
      
      # curve(dgamma(x, shape=alpha, rate=lambda), from=0, to=15, n=500, type="l", xlab = "Wind Speed", ylab = "")
      # abline(v = observed)
      # abline(v = percent[1], col = 'blue')
      # abline(v = percent[2], col = 'blue')
    
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
  score = count/(3 * diff)
  scores = c(scores, score)
  print(score)
}
# plot(diffs, scores, type="o")
# abline(h = upperBound - lowerBound, col = 'blue')