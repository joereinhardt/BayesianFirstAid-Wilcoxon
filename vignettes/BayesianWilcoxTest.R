## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
pairDiff_means <- function(N) {
  #create x and y so that x > y
  x <- seq(1.1, 2.1, by = (1/(N - 1)))
  y <- seq(0, 1, by = (1/(N - 1)))
  n <- length(x) + length(y)
  ranks <- rank(c(x, y)) # Replace by ranks
  #Tranformed Ranks
  seqQ <- seq(1, 2*n - 1, by = 2)/(2*n)
  zRanks <- qnorm(seqQ[ranks])
  zRanksX <- zRanks[1:length(x)]
  zRanksY <- zRanks[-(1:length(x))]
  mean_diff <- mean(zRanksX - zRanksY)
  return(mean_diff)
}
means <- c()
for (i in 2:500) {
  means <- c(means, pairDiff_means(i))
}
plot(means, type = "l", ylab = "mu_max",  xlab = "N")

