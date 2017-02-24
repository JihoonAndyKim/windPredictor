windData = testWind
windDir = testDirWind

end = dim(windData)[1]
startWindow = 1
endWindow = 40
windowLength = endWindow

timeSeries <- vector()
while(endWindow < end) {
  
  m <- windData[startWindow:endWindow, 1:11, 1]
  x <- vector()
  for(i in seq(nrow(m))) {
   x <- c(x, mean(m[i, 1:11])) 
  }
  y <- windData[startWindow:endWindow, 12, 1]
  startWindow = startWindow + 1
  endWindow = endWindow + 1
  df = data.frame(x, y)
  
  timeSeries <- c(timeSeries, cor(df, method = "pearson")[2])
}

plot(1:length(timeSeries), timeSeries, type = 'l', col = "red")
abline(h = 0)
abline(v = seq(1,7000, by = 365))
#plot(1:1000, timeSeries[1:1000])
highCorrelation <- which(timeSeries[1:1000] >= 0.6)

#timeSlice = as.integer(highCorrelation[1])
n = windDir[818:858, 1:11,1]
avgDir <- vector()
for(i in seq(nrow(n))) {
  mean = mean(n[i, 1:11])
  if(mean < 0) {
    mean = mean + 360
  }
  avgDir <- c(avgDir, mean) 
}

#direction = c("N","NNE","NE","ENE","E","ESE", "SE", "SSE","S","SSW","SW","WSW","W","WNW","NW","NNW")
#avgDir <- lapply(avgDir, function(p){ direction[as.integer((p/22.5)+.5) %% 16]})

m <- windData[818:858,,1]
avgWind <- vector()
obsWind <- m[, 12]
for(i in seq(nrow(m))) {
  avgWind <- c(avgWind, mean(m[i, 1:11])) 
}

plot(avgWind, type = "o", col = "red")
#plot the Observation data
lines(obsWind, type = "o")

#see if wind direction of influence (Or is it seasonal?) Try to understand this
#spec.pgram (dft with squared coefficients)
#45 degree bins instead. Check for correlation within bins
#rose diagram --> correlation is a function of direction.Between obs and forecast
#Plot season
#do 2 and 3 day ahead as