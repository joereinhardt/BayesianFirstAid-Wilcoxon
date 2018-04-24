#' Bayesian First Aid Alternative to the Wilcoxon Rank Test
#'
#' Implements a Bayesian alternative to the \code{\link{wilcox.test}} function.
#' Written as an addition to and based on the Bayesian First Aid package by
#' Rasmus Bååth (2014). The method follows
#' Gelman et al. (2014), and
#' fits a normal model using JAGS after an inverse-normal rank-transformation.
#'
#' The test performs an inverse-normal rank transformation on the data before
#' the test models are run. This transformation is based on Gelman et al. (2014)
#' and transforms the ranks of the data to quantiles of a standard Gaussian. See
#' the package vignette for more details.
#'
#' If \code{paired = FALSE}, an independent (two) sample Wilcoxon test is
#' performed.
#' For the two-sample wilcox test, the posterior distribution is obtained from
#' the following model:
#'
#' \figure{twoSampleWilcoxDiagram.svg}{options: height=400}
#'
#' If \code{paired = TRUE}, a paired sample test is run, with the following
#' model structure:
#'
#' \figure{oneSampleWilcoxDiagram.svg}{options: height=275}
#'
#' If \code{mu} is given, the hypothesis that the inverse-normal rank
#' transformed data are larger/smaller than \code{mu} is tested.
#'
#' @param x numeric vector of data values
#' @param y numeric vector of data values to be compared to x. Unlike for the
#'   \code{\link{wilcox.test}} function, this argument is required.
#' @param cred.mass the amount of probability mass that will be contained in
#'   reported credible intervals. This argument fills a similar role as
#'   conf.level in \code{\link{wilcox.test}}
#' @param mu number specifying an optional paramter to form the null hypothesis.
#'   See 'Details'.
#' @param paired a logical indicating whether you want a paired test.
#' @param n.iter The number of iterations to run the MCMC sampling.
#' @param alternative ignored, only retained in order to maintain compability
#'   with \code{\link{wilcox.test}}
#' @param exact ignored, only retained in order to maintain compability with
#'   \code{\link{wilcox.test}}
#' @param correct ignored, only retained in order to maintain compability with
#'   \code{\link{wilcox.test}}
#' @param conf.int ignored, only retained in order to maintain compability with
#'   \code{\link{wilcox.test}}
#' @param conf.level identical to cred.mass, ignored, only retained in order to
#'   maintain compability with \code{\link{wilcox.test}}
#' @param progress.bar The type of progress bar. Possible values are "text",
#'   "gui", and "none".
#' @param ... further arguments to be passed to or from methods.
#'
#' @return A list of class \code{bayes_wilcox_test}. It can be further inspected
#'   using the functions \code{print}, \code{summary}, \code{plot},
#'   \code{\link{diagnostics.bayes_paired_wilcox_test}} and
#'   \code{\link{model.code.bayes_paired_wilcox_test}}, or
#'   \code{\link{diagnostics.bayes_two_sample_wilcox_test}} and
#'   \code{\link{model.code.bayes_two_sample_wilcox_test}}, respectively.
#'
#' @examples
#' # Using examples copied from the wilcox.test help file:
#'
#' # One-sample (Paired sample) test.
#' # Hollander & Wolfe (1973), 29f.
#' # Hamilton depression scale factor measurements in 9 patients with
#' # mixed anxiety and depression, taken at the first (x) and second
#' # (y) visit after initiation of a therapy (administration of a
#' # tranquilizer).
#' x <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
#' y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
#' # Performing a paired-sample Bayesian Wilcox test. The function can be used
#' # in exactly the same way as wilcox.test:
#' pairedBayesWilcox <- bayes.wilcox.test(x, y, paired = TRUE,
#'       alternative = "greater")
#' # Note that the argument "alternative" is ignored by bayes.wilcox.test
#' # Test summary
#' summary(pairedBayesWilcox)
#' # Visual Inspection
#' plot(pairedBayesWilcox)
#' # Diagnostics for MCMC
#' diagnostics.bayes_paired_wilcox_test(pairedBayesWilcox)
#' # Print out the model code, which can be modified
#' model.code.bayes_paired_wilcox_test(pairedBayesWilcox)
#' # Classical Wilcox Test for comparison
#' wilcox.test(x, y, paired = TRUE, alternative = "greater")
#'
#'
#' # Two-sample (independent) test.
#' # Hollander & Wolfe (1973), 69f.
#' # Permeability constants of the human chorioamnion (a placental
#' # membrane) at term (x) and between 12 to 26 weeks gestational
#' # age (y).  The alternative of interest is greater permeability
#' # of the human chorioamnion for the term pregnancy.
#' x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
#' y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
#' indBayesWilcox <- bayes.wilcox.test(x, y, alternative = "greater")
#'                               # Note that the argument "alternative" is
#'                               # ignored by bayes.wilcox.test
#' # Test summary
#' summary(indBayesWilcox)
#' # Visual Inspection
#' plot(indBayesWilcox)
#' # Diagnostics for MCMC
#' diagnostics.bayes_two_sample_wilcox_test(indBayesWilcox)
#' # Print out the model code, which can be modified
#' model.code.bayes_two_sample_wilcox_test(indBayesWilcox)
#' # Classical Wilcox Test for comparison
#' wilcox.test(x, y, alternative = "greater")
#'
#' @references Rasmus Bååth (2014), “Bayesian First Aid: A Package That
#'   Implements Bayesian Alternatives to the Classical *.Test Functions in R.”
#'   In \emph{UseR! 2014 - the International R User Conference.}
#' @references Andrew Gelman, John B Carlin, Hal S Stern,
#'  David B Dunson, Aki Vehtari and Donald B Rubin (2014),
#'  \emph{Bayesian Data Analysis}. 3rd ed. CRC press Boca Raton, FL. Page 97.
#' @references Myles Hollander and Douglas A. Wolfe (1973), \emph{Nonparametric
#'   Statistical Methods.} New York: John Wiley & Sons. Pages 27–33
#'   (one-sample), 68–75 (two-sample). Or second edition (1999).

#'
#' @import coda
#' @import rjags
#' @import stats
#' @import graphics
#' @import grDevices
#' @import MASS
#'
#' @export
#' @rdname bayes.wilcox.test

bayes.wilcox.test <- function(x, ...){
  UseMethod("bayes.wilcox.test")
}


#' @export
#' @rdname bayes.wilcox.test
bayes.wilcox.test.default <- function(x, y, cred.mass = 0.95,
                                      mu = 0,
                                      paired = FALSE,
                                      n.iter = 30000,
                                      alternative = NULL,
                                      exact = NULL,
                                      correct = NULL,
                                      conf.int = NULL,
                                      conf.level,
                                      progress.bar = "text", ...) {
  ### Original (but slighly modified) code from t.test.default ###

  #If statements from BFA
  if (!missing(conf.level)) {
    cred.mass <- conf.level
  }

  if (!missing(alternative)) {
    warning("The argument 'alternative' is ignored by bayes.wilcox.test")
  }

  if (!missing(exact)) {
    warning("The argument 'exact' is ignored by bayes.wilcox.test")
  }

  if (!missing(correct)) {
    warning("The argument 'correct' is ignored by bayes.wilcox.test")
  }

  if (!missing(conf.int)) {
    warning("The argument 'conf.int' is ignored by bayes.wilcox.test")
  }


  ### Original (but slighly modified) code from t.test.default ###

  if (!missing(mu) && (length(mu) != 1 || is.na(mu)))
      stop("'mu' must be a single number")
  if (!missing(cred.mass) && (length(cred.mass) != 1 || !is.finite(cred.mass) ||
                                cred.mass < 0 || cred.mass > 1))
      stop("'cred.mass' or 'conf.level' must be a single number between 0 and 1")


  # removing incomplete cases and preparing the data vectors (x & y)
  x_name <- deparse(substitute(x))
  y_name <- deparse(substitute(y))
  data_name <- paste(x_name, "and", y_name)

  if (paired)
    xok <- yok <- complete.cases(x, y)
  else {
    yok <- !is.na(y)
    xok <- !is.na(x)
  }
  y <- y[yok]
  x <- x[xok]
  # Checking that there is enough data
  nx <- length(x)
  ny <- length(y)
  n <- nx + ny
  if (nx < 2)
    stop("not enough 'x' observations")
  if (ny < 2)
    stop("not enough 'y' observations")

  #Transform Data to Inverse Normal Ranks
  ranks <- rank(c(x, y))
  ranksX <- ranks[1:nx]
  ranksY <- ranks[-(1:nx)]
  seqQ <- seq(1, 2*n - 1, by = 2)/(2*n)
  qRanks <- seqQ[ranks]
  zRanks <- qnorm(qRanks)
  zRanksX <- zRanks[1:nx]
  zRanksY <- zRanks[-(1:nx)]

  ### Running Model using JAGS. Code adapted from BFA
  #Run model
  if (paired == TRUE) {
  mcmc_samples <- jags_paired_wilcox_test(zRanksX, zRanksY,
                                          n.chains = 3,
                                          n.iter = ceiling(n.iter / 3),
                                          progress.bar = progress.bar)
  stats <- mcmc_stats(mcmc_samples, cred_mass = cred.mass,
                      comp_val = mu) #mcmc_stats function from BFA's utilities
  #quick fix to get around utility function not assigning names if only one var
  rownames(stats) <- "mu_diff"
  bfa_object <- list(x = x, y = y, pair_diff = x - y, comp = mu,
                     cred_mass = cred.mass, x_name = x_name, y_name = y_name,
                     data_name = data_name, x_data_expr = x_name,
                     y_data_expr = y_name, mcmc_samples = mcmc_samples,
                     stats = stats)
  class(bfa_object) <- c("bayes_paired_wilcox_test", "bayesian_first_aid")

  } else {
  mcmc_samples <- jags_two_sample_wilcox_test(zRanksX, zRanksY, n.chains = 3,
                                         n.iter = ceiling(n.iter / 3),
                                         progress.bar = progress.bar)
  stats <- mcmc_stats(mcmc_samples, cred_mass = cred.mass, comp_val = mu)
  #mcmc_stats function from BFA's utilities
  bfa_object <- list(x = x, y = y, comp = mu, cred_mass = cred.mass,
                     x_name = x_name, y_name = y_name, data_name = data_name,
                     x_data_expr = x_name, y_data_expr = y_name,
                     mcmc_samples = mcmc_samples, stats = stats)
  class(bfa_object) <- c("bayes_two_sample_wilcox_test", "bayesian_first_aid")

  }
  bfa_object
}

#Functions for running JAGS are outsourced to a separate function.

#Paired-Sample model
####################
#JAGS model string
paired_samples_wilcox_model_string <- "model {
for (i in 1:length(pair_diff)) {
  pair_diff[i] ~ dnorm(mu_diff, sigma_diff)
}
mu_diff ~ dunif(-1.6, 1.6)
sigma_diff ~ dunif(0, sigma_high)
}"



#Function for paired JAGS model
jags_paired_wilcox_test <- function(x, y, n.adapt = 1000,
                                    n.chains = 3, n.update = 500,
                                    n.iter = 5000, thin = 1,
                                    progress.bar = "text") {
  pair_diff <- x - y
  n <- length(x) + length(y)
  #Standard JAGS structure.
  data_list <- list(pair_diff = pair_diff,
                    sigma_high = -2*qnorm(1/(2*n)))
  inits_list <- list(mu_diff = 0, sigma_diff = 1)
  params <- c("mu_diff")
  #Using run_jags function from BFA utilities
  mcmc_samples <- run_jags(paired_samples_wilcox_model_string,
                           data = data_list,
                           inits = inits_list,
                           params = params,
                           n.chains = n.chains,
                           n.adapt = n.adapt,
                           n.update = n.update,
                           n.iter = n.iter,
                           thin = thin,
                           progress.bar = progress.bar)
  mcmc_samples

}

#Two-sample model
#################

#JAGS model string
two_sample_wilcox_model_string <- "model {
for(i in 1:length(x)) {
  x[i] ~ dnorm( mu_x , sigma_x )
}
for(i in 1:length(y)) {
  y[i] ~ dnorm( mu_y , sigma_y )
}
mu_diff <- mu_x - mu_y

# The priors
mu_x ~ dunif(mu_low, mu_high)
sigma_x ~ dunif(0, sigma_high)
mu_y ~ dunif(mu_low, mu_high)
sigma_y ~ dunif(0, sigma_high)
}"

#Function for two sample test, using run_jags function from BFA utilities
jags_two_sample_wilcox_test <- function(x, y,
                                        n.adapt= 1000,
                                        n.chains = 3, n.update = 500,
                                        n.iter = 10000, thin = 1,
                                        progress.bar = "text") {
  n <- length(x) + length(y)
  data_list <- list(x = x,
                    y = y,
                    mu_low = qnorm(1/(2*n)),
                    mu_high = qnorm((2*n - 1)/(2*n)),
                    sigma_high = -qnorm(1/(2*n)))

  inits_list <- list( #Assume parameters under equal distribution for x, y
    mu_x = 0,
    mu_y = 0,
    sigma_x = 1,
    sigma_y = 1)

  params <- c("mu_x", "mu_y", "mu_diff")
  mcmc_samples <- run_jags(two_sample_wilcox_model_string,
                           data = data_list,
                           inits = inits_list,
                           params = params,
                           n.chains = n.chains, n.adapt = n.adapt,
                           n.update = n.update, n.iter = n.iter,
                           thin = thin, progress.bar = progress.bar)
  mcmc_samples
}

############################################
#S3 functions. All S3 functions are adapted
#from BFA's bayes.t.test or bayes.binom.test
#S3 functions

############################################
### Paired Sample Wilcox Test S3 Methods ###
############################################

#' @export
# Concise Information printout
print.bayes_paired_wilcox_test <- function(x, ...) {
  s <- format_stats(x$stats)

  cat("\n")
  cat("\tBayesian First Aid Wilcoxon test\n")
  cat("\n")

  cat("data: ", x$x_name, " (n = ", length(x$x) ,") and ",
      x$y_name, " (n = ", length(x$y) ,")\n", sep = "")
  cat("\n")

  cat("  Estimates [", s[1, "HDI%"] ,"% credible interval]\n", sep = "")
  cat("mean paired difference: ", s["mu_diff", "median"],
      " [", s["mu_diff", "HDIlo"], ", ",
      s["mu_diff", "HDIup"] , "]\n", sep = "")
  cat("\n")
  cat("The mean difference is more than",
      s["mu_diff","comp"] , "by a probability of", s["mu_diff","%>comp"], "\n")
  cat("and less than", s["mu_diff", "comp"],
      "by a probability of", s["mu_diff", "%<comp"], "\n")
  cat("\n")
  invisible(NULL)
}

#' @method summary bayes_paired_wilcox_test
#' @export
# More detailed summary printout
summary.bayes_paired_wilcox_test <- function(object, ...) {
  s <- round(object$stats, 3)

  cat("  Data\n")
  cat(object$x_name, ", n = ", length(object$x), "\n", sep = "")
  cat(object$y_name, ", n = ", length(object$y), "\n", sep = "")
  cat("\n")

  #print_bayes_two_sample_t_test_params(object) #Replace by shorter:
  cat("  Model parameters and generated quantities\n")
  cat("mu_diff: the difference in means of the inverse-normal transformed ranks
      of", object$x_name, "and",
      object$y_name, "\n")
  cat("\n")

  cat("  Measures\n" )
  print(s[, c("mean", "HDIlo", "HDIup", "%<comp", "%>comp")])
  cat("\n")
  cat("'HDIlo' and 'HDIup' are the limits of a ", s[1, "HDI%"] ,
      "% HDI credible interval.\n", sep = "")
  cat("'%<comp' and '%>comp' are the probabilities of the respective parameter",
      "being\n")
  cat("smaller or larger than ", s[1, "comp"] ,".\n", sep = "")

  cat("\n")
  cat("  Quantiles\n" )
  print(s[, c("q2.5%", "q25%", "median","q75%", "q97.5%")] )
  invisible(object$stats)
}

#adapted from plot.bayes_binom_test. Plots the posterior MCMC samples
#' @export
plot.bayes_paired_wilcox_test <- function(x, ...) {
  old_par <- par(mar = c(3.5,3.5,2.5,0.5),
                 mgp = c(2.25,0.7,0),
                 mfcol = c(1,1))
  stats <- x$stats
  mcmc_samples <- x$mcmc_samples
  samples_mat <- as.matrix(mcmc_samples)
  mu = 0
  sample_mat <- as.matrix(x$mcmc_samples)
  xlim = range(c(sample_mat[,"var1"], x$comp))
  plotPost(sample_mat[,"var1"], cred_mass = x$cred_mass,
           comp_val = x$comp, xlim = xlim, cex = 1, cex.lab = 1.5,
           main = "Difference in means of the quantile-rank transformed data",
           xlab = bquote(mu[diff]), show_median = TRUE)
  par(old_par)
  invisible(NULL)
}

#Adapted from model.code.bayes_binom_test
#' Prints code that replicates the model you just ran.
#'
#' This is good if you better want to understand how the model is
#' implemented or if you want to run a modified version of the code.
#'
#' @export
#' @param fit The output from a paired sample Bayesian Wilcox test model.


# Note: function not recognized as S3, needs to be called by
# model.code.bayes_paired_wilcox_test().
# #' method model.code bayes_paired_wilcox_test
# does not fix the issue. Same problem for:
# model.code.bayes_two_sample_wilcox_test,
# diagnostics.bayes_paired__wilcox_test,
# diagnostics.bayes_two_sample_wilcox_test


model.code.bayes_paired_wilcox_test <- function(fit) {
  cat("### Model code for the Bayesian First Aid alternative to the Wilcox",
      "paired sample test ###\n\n")
  cat("require(rjags)\n\n")

  cat("# Setting up the data\n")
  cat("x <-", fit$x_name, "\n")
  cat("y <-", fit$y_name, "\n")
  cat("\n")
  pretty_print_function_body(paired_samples_wilcox_model_code)
  invisible(NULL)
}

# Not to be run, just to be printed - adapted from two_sample_t_test_model_code
paired_samples_wilcox_model_code <- function(pair_diff) {
  # The model string written in the JAGS language
  BayesianFirstAid::replace_this_with_model_string

  # Initializing parameters to sensible starting values helps the convergence
  # of the MCMC sampling. Here using robust estimates of the mean (trimmed)
  # and standard deviation (MAD).
  inits_list <- list(mu_diff = mean(pair_diff, trim = 0.2))

  n <- length(pair_diff)
  #Standard JAGS structure.
  data_list <- list(pair_diff = pair_diff,
                    sigma_high = -2*qnorm(1/(2*n)))
  inits_list <- list(mu_diff = 0, sigma_diff = 1)


  # The parameters to monitor.
  params <- c("mu_diff")

  # Running the model
  model <- jags.model(textConnection(model_string), data = data_list,
                      inits = inits_list, n.chains = 3, n.adapt = 1000)
  update(model, 500) # Burning some samples to the MCMC gods....
  samples <- coda.samples(model, params, n.iter = 10000)

  # Inspecting the posterior
  plot(samples)
  summary(samples)
}

paired_samples_wilcox_model_code <- inject_model_string(
  paired_samples_wilcox_model_code,
  paired_samples_wilcox_model_string)


#adapted from diagnostics.bayes_two_sample_t_test
#' Plots and prints diagnostics regarding the convergence of the model.
#'
#' @param x The output from a paired sample Bayesian Wilcox test model.
#' @export

diagnostics.bayes_paired_wilcox_test <- function(x) {
  print_mcmc_info(x$mcmc_samples)
  cat("\n")
  print_diagnostics_measures(round(x$stats, 3))
  cat("\n")
  cat("Model parameters and generated quantities\n")
  cat("mu_diff: the difference in means of the quantile-rank transformed data",
      x$x_name, "and",
      x$y_name, "\n")

  cat("\n")

  old_par <- par(mar = c(3.5,2.5,2.5,0.5) , mgp = c(2.25,0.7,0))
  plot(x$mcmc_samples)
  par(old_par)
  invisible(NULL)
}

### End of Paired Sample Wilcox Test S3 Functions ###

#########################################
### Two Sample Wilcox Test S3 Methods ###
#########################################

#Concise printout.
#' @export
print.bayes_two_sample_wilcox_test <- function(x, ...) {
  s <- format_stats(x$stats)

  cat("\n")
  cat("\tBayesian First Aid Wilcoxon test\n")
  cat("\n")

  cat("data: ", x$x_name, " (n = ", length(x$x) ,") and ",
      x$y_name, " (n = ", length(x$y) ,")\n", sep = "")
  cat("\n")

  cat("  Estimates [", s[1, "HDI%"] ,"% credible interval]\n", sep = "")
  cat("mean of ",  x$x_name, ": ", s["mu_x", "median"],
      " [", s["mu_x", "HDIlo"], ", ", s["mu_x", "HDIup"] , "]\n",sep = "")
  cat("mean of ",  x$y_name, ": ", s["mu_y", "median"],
      " [", s["mu_y", "HDIlo"], ", ", s["mu_y", "HDIup"] , "]\n",sep = "")
  cat("difference of the means: ", s["mu_diff", "median"],
      " [", s["mu_diff", "HDIlo"], ", ", s["mu_diff", "HDIup"] , "]\n",sep = "")

  cat("\n")
  cat("The mean difference is more than",
      s["mu_diff","comp"] , "by a probability of", s["mu_diff","%>comp"], "\n")
  cat("and less than", s["mu_diff", "comp"],
      "by a probability of", s["mu_diff", "%<comp"], "\n")
  cat("\n")
  invisible(NULL)
}

#More detailed summary printout.
#' @method summary bayes_two_sample_wilcox_test
#' @export
summary.bayes_two_sample_wilcox_test <- function(object, ...) {
  s <- round(object$stats, 3)

  cat("  Data\n")
  cat(object$x_name, ", n = ", length(object$x), "\n", sep = "")
  cat(object$y_name, ", n = ", length(object$y), "\n", sep = "")
  cat("\n")

  #print_bayes_two_sample_t_test_params(object) #Replace by shorter:
  cat("  Model parameters and generated quantities\n")
  cat("mu_x: the mean of the quantile-rank transformed", object$x_name, "\n")
  cat("mu_diff: the mean of the quantile-rank transformed ", object$x_name, "and",
      object$y_name, "\n")
  cat("mu_y: the mean of the quantile-rank transformed", object$y_name, "\n")
  cat("\n")

  cat("  Measures\n" )
  print(s[, c("mean", "HDIlo", "HDIup", "%<comp", "%>comp")])
  cat("\n")
  cat("'HDIlo' and 'HDIup' are the limits of a ", s[1, "HDI%"] ,
      "% HDI credible interval.\n", sep = "")
  cat("'%<comp' and '%>comp' are the probabilities of the respective parameter",
      "being\n")
  cat("smaller or larger than ", s[1, "comp"] ,".\n", sep = "")

  cat("\n")
  cat("  Quantiles\n" )
  print(s[, c("q2.5%", "q25%", "median","q75%", "q97.5%")] )
  invisible(object$stats)
}

#Plots the posterior MCMC samples
#' @export
plot.bayes_two_sample_wilcox_test <- function(x, ...) {
  stats <- x$stats
  mcmc_samples <- x$mcmc_samples
  samples_mat <- as.matrix(mcmc_samples)
  mu_x = samples_mat[,"mu_x"]
  mu_y = samples_mat[,"mu_y"]
  mu_diff = samples_mat[,"mu_diff"]

  old_par <- par(no.readonly = TRUE)
  layout(matrix(c(3, 1,
                  3, 2) , nrow = 2, byrow = FALSE))
  par(mar = c(3.5,3.5,2.5,0.51) , mgp = c(2.25,0.7,0) )

  # Plot posterior distribution of parameters of the quantile-rank transformed
  # mu_x, mu_y, and their difference:
  xlim = range(c(mu_x , mu_y))
  plotPost(mu_x,  xlim = xlim, cex.lab = 1.75, cred_mass = x$cred_mass,
           show_median = TRUE, xlab = bquote(mu[x]),
           main = paste(x$x_name,"Mean"), col = "skyblue")
  plotPost(mu_y,  xlim = xlim, cex.lab = 1.75,  cred_mass = x$cred_mass,
           show_median = TRUE, xlab = bquote(mu[y]),
           main = paste(x$y_name,"Mean") , col = "skyblue")
  plotPost(mu_diff, comp_val = x$comp, cred_mass = x$cred_mass,
            xlab = bquote(mu[x] - mu[y]) , cex.lab = 1.75 , show_median = TRUE,
            main = "Difference of Means" , col = "skyblue")

  par(old_par)
  invisible(NULL)
  }

#' Prints code that replicates the model you just ran.
#'
#' This is good if you better want to understand how the model is
#' implemented or if you want to run a modified version of the code.
#'
#' @param fit The output from a two sample Bayesian Wilcox test model.
#' @export
model.code.bayes_two_sample_wilcox_test <- function(fit) {
  cat("### Model code for the Bayesian First Aid alternative to the Wilcox",
      "two-sample test ###\n\n")
  cat("require(rjags)\n\n")

  cat("# Setting up the data\n")
  cat("x <-", fit$x_name, "\n")
  cat("y <-", fit$y_name, "\n")
  cat("\n")
  pretty_print_function_body(two_sample_wilcox_model_code)
  invisible(NULL)
}

# Not to be run, just to be printed - adapted from two_sample_t_test_model_code
two_sample_wilcox_model_code <- function(x, y) {
  # The model string written in the JAGS language
  BayesianFirstAid::replace_this_with_model_string

  # Initializing parameters to sensible starting values helps the convergence
  # of the MCMC sampling. Here using parameters for both samples having the
  # same distribution
  inits_list <- list(
    mu_x = 0,
    mu_y = 0,
    sigma_x = 1,
    sigma_y = 1)

  #Setting parameters for the priors
  n <- length(x) + length(y)
  data_list <- list(x = x,
                    y = y,
                    mu_low = qnorm(1/(2*n)),
                    mu_high = qnorm((2*n - 1)/(2*n)),
                    sigma_high = -qnorm(1/(2*n)))

  # The parameters to monitor.
  params <- c("mu_x", "mu_y", "mu_diff")

  # Running the model
  model <- jags.model(textConnection(model_string), data = data_list,
                      inits = inits_list, n.chains = 3, n.adapt = 1000)
  update(model, 500) # Burning some samples to the MCMC gods....
  samples <- coda.samples(model, params, n.iter = 5000)

  # Inspecting the posterior
  plot(samples)
  summary(samples)
}

two_sample_wilcox_model_code <- inject_model_string(
  two_sample_wilcox_model_code,
  two_sample_wilcox_model_string)
#adapted from diagnostics.bayes_two_sample_t_test
#' Plots and prints diagnostics regarding the convergence of the model.
#'
#' @param x The output from a two sample Bayesian Wilcox test model.
#' @export
# Note: function not recognized as S3, needs to be called by
# model.code.bayes_wilcox_test()

diagnostics.bayes_two_sample_wilcox_test <- function(x) {
  print_mcmc_info(x$mcmc_samples)
  cat("\n")
  print_diagnostics_measures(round(x$stats, 3))
  cat("\n")
  cat("  Model parameters and generated quantities\n")
  cat("mu_x: the mean of the quantile-rank transformed", x$x_name, "\n")
  cat("mu_diff: the difference in means of the quantile-rank
      transformed", x$x_name, "and",
      x$y_name, "\n")
  cat("mu_y: the mean of the quantile-rank transformed", x$y_name, "\n")
  cat("\n")

  old_par <- par( mar = c(3.5,2.5,2.5,0.5) , mgp = c(2.25,0.7,0) )
  plot(x$mcmc_samples)
  par(old_par)
  invisible(NULL)
}


### End of Two Sample Wilcox Test S3 Functions ###
