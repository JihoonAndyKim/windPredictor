
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



