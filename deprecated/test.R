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

obs <- formatDate(read.csv("WindData.csv"))