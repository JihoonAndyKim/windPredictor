library(mvtnorm)

genBivariate = abs(rmvnorm(n = 1000, mean = c(0, 0), sigma = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)))

plot(dgamma(abs(genBivariate[,1]), shape = 3, scale = 1), dgamma(abs(genBivariate[,2]), shape = 3, scale = 1))