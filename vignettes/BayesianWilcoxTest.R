## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE----------------------------------------------------
require(BayesianFirstAid)

## ------------------------------------------------------------------------
#Setting up example data
data <- c(68, 24)
#Classical test
binom.test(data, p = 3/4, alternative = "l")

## ---- fig.width = 5, fig.height = 5--------------------------------------
#Bayesian Alternative
bayesBinom <- bayes.binom.test(data, p = 3/4)
print(bayesBinom)
plot(bayesBinom)

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

## ----graphPaired, echo=FALSE, out.width = '55%'--------------------------
knitr::include_graphics("oneSampleWilcoxDiagram.svg")

## ----graphTwo, echo=FALSE, out.width = '75%'-----------------------------
knitr::include_graphics("twoSampleWilcoxDiagram.svg")

## ----eval = FALSE--------------------------------------------------------
#  require(devtools)
#  install_github("joereinhardt/BayesianFirstAid-Wilcoxon")

## ------------------------------------------------------------------------
#Copied from the wilcox.test help file:
# Hamilton depression scale factor measurements in 9 patients with
#  mixed anxiety and depression, taken at the first (x) and second
#  (y) visit after initiation of a therapy (administration of a
#  tranquilizer).
x <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
wilcox.test(x, y, paired = TRUE, alternative = "greater")

## ---- fig.width = 5------------------------------------------------------
require(bayesWilcoxTest)
hamiltonBayesWilcox <- bayes.wilcox.test(x, y, paired = TRUE,
                                         alternative = "greater")

#Print out concise information on the test
print(hamiltonBayesWilcox)

#Show a more detailed summary
summary(hamiltonBayesWilcox)

#Visual Inspection of the posterior
plot(hamiltonBayesWilcox)

#MCMC diagnostics
diagnostics.bayes_paired_wilcox_test(hamiltonBayesWilcox)

#Obtain model code in order to modify as needed
model.code.bayes_paired_wilcox_test(hamiltonBayesWilcox)

## ------------------------------------------------------------------------
#Copied from the wilcox.test help file:
# Permeability constants of the human chorioamnion (a placental
#  membrane) at term (x) and between 12 to 26 weeks gestational
#  age (y).  The alternative of interest is greater permeability
#  of the human chorioamnion for the term pregnancy.
x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
wilcox.test(x, y, alternative = "g")

## ---- fig.width = 5, fig.height = 5--------------------------------------
membraneBayesWilcox <- bayes.wilcox.test(x, y, alternative = "g")

#Print out concise information on the test
print(membraneBayesWilcox)

#Show a more detailed summary
summary(membraneBayesWilcox)

#Visual Inspection of the posterior
plot(membraneBayesWilcox)

#MCMC diagnostics
diagnostics.bayes_two_sample_wilcox_test(membraneBayesWilcox)

#Obtain model code in order to modify as needed
model.code.bayes_two_sample_wilcox_test(membraneBayesWilcox)

