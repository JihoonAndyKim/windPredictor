# formatDate <- function(obs) {
#   obs <- obs[obs$Hour == 0,]
#   obs$Date <- gsub("-", "", paste(as.Date(obs$Date, "%m/%d/%Y"), "00", sep= ""))
#   obs$Wind.Speed <- as.numeric(as.character(obs$Wind.Speed))
#   return(obs)
# }
# 
# obs <- formatDate(read.csv("WindData.csv"))


vec <- vector()

for(i in 7070 : 7101) {
  vec <- c(vec, mean(testWind[i, 1:11, 1]))
}

# Plot the Ensemble member forecast
plot(vec, type = "o", col = "red")
#plot the Observation data
lines(obs[7070:7101,3], type = "o")

# plot(vec, obs[6101:7101,3])
# print(cor(vec, obs[6101:7101, 3], use = "pairwise"))

# plot(testWind[1:50, 1, 1], type = "o")
# for(i in 2:11) {
#   lines(testWind[1:50, i, 1], type = "o")
# }

#ccf

