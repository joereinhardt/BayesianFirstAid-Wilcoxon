#' Bayesian First Aid Alternative to the t-test
#'
#' Function implements a bayesian alternative to the wilcoxon.test() function.
#' Written as an addition to the Bayesian First Aid (BFA) package by
#' Rasmus Bååth (https://github.com/rasmusab). BFA's structure is largely
#' followed. S3 functions to be implemented as in BFA:
#' print, plot, summary, model.code, and diagnostics.
#' Method follows http://andrewgelman.com/2015/07/13/dont-do-the-wilcoxon/, and
#' fits a normal model using JAGS after rank-transformation.
#'
#' @param x numeric vector of data values
#' @param y numeric vector of data values to be compared to x
#' @param cred.mass he amount of probability mass that will be contained in
#'      reported credible intervals. This argument fills a similar role as
#'      conf.level in \code{\link{wilcox.test}}
#' @param mu a number specifying an optional parameter used to form the null
#'     hypothesis.
#' @param paired a logical indicating whether you want a paired test.
#' @param n.iter The number of iterations to run the MCMC sampling.
#' @param alternative ignored, only retained in order to maintain compability
#'     with \code{\link{wilcox.test}}
#' @param conf.level identical to cred.mass,
#'     ignored, only retained in order to maintain compability with
#'     \code{\link{wilcox.test}}
#' @param progress.bar The type of progress bar. Possible values are "text",
#'  "gui", and "none".
#' @param ... further arguments to be passed to or from methods.
#'
#' @return A list of class \code{bayes_wilcox_test}. It can be further inspected
#'  using the functions \code{summary}, \code{plot}, \code{diagnostics}
#'  and \code{model.code}.
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
  UseMethod("bayes.t.test")
}

#Note: needs to take same input as wilcox.test()
#' @export
#' @rdname bayes.wilcox.test
bayes.wilcox.test.default <- function(x, y, cred.mass = 0.95,
                                      mu = 0,
                                      paired = FALSE,
                                      n.iter = 30000,
                                      alternative = NULL,
                                      conf.level,
                                      progress.bar = "text", ...) {
  ### Original (but slighly modified) code from t.test.default ###

  #If statements from BFA
  if (!missing(conf.level)) {
    cred.mass <- conf.level
  }

  if (!missing(alternative)) {
    warning("The argument 'alternative' is ignored by bayes.binom.test")
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
  if (nx < 2)
    stop("not enough 'x' observations")
  if (ny < 2)
    stop("not enough 'y' observations")

  ### Running Model. Code adapted from BFA

  if (paired) {
  mcmc_samples <- jags_paired_wilcox_test(x, y,
                                          n.chains = 3,
                                          n.iter = ceiling(n.iter / 3),
                                          progress.bar = progress.bar)
  stats <- mcmc_stats(mcmc_samples, cred_mass = cred.mass,
                      comp_val = mu)
  bfa_object <- list(x = x, y = y, pair_diff = x - y, comp = mu,
                     cred_mass = cred.mass, x_name = x_name, y_name = y_name,
                     data_name = data_name, x_data_expr = x_name,
                     y_data_expr = y_name, mcmc_samples = mcmc_samples,
                     stats = stats)
  class(bfa_object) <- c("bayes_paired_wilcox_test", "bayesian_first_aid")

  } else {
  NULL
  }
  bfa_object
}

#JAGS model string
paired_samples_wilcox_model_string <- "model {
  for (i in 1:length(pair_diff)) {
    pair_diff[i] ~ dnorm(mu_diff, 1)
  }

  mu_diff ~ dunif(-1.6, 1.6)
}"

#Figure out how to include comp.mu!

#Function for JAGS model
jags_paired_wilcox_test <- function(x, y, comp.mu = 0, n.adapt = 500,
                                    n.chains = 3, n.update = 100,
                                    n.iter = 5000, thin = 1,
                                    progress.bar = "text") {
  pair_diff <- x - y
  data_list <- list(
    pair_diff = pair_diff,
    comp.mu = comp.mu
  )

  inits_list <- list(mu_diff = mean(pair_diff, trim = 0.2))
  params <- c("mu_diff")
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

### Paired Sample Wilcox Test S3 Methods ###


#' @export
print.bayes_paired_wilcox_test <- function(x, ...) {
  s <- format_stats(x$stats)

  cat("\n")
  cat("\tBayesian First Aid Wilcoxon test\n")
  cat("\n")
  cat("data: ", x$x_name, " (n = ", length(x$x) ,") and ",
      x$y_name, " (n = ", length(x$y) ,")\n", sep = "")
  cat("\n")
  cat("  Estimates [", s[1, "HDI%"] ,"% credible interval]\n", sep = "")
  cat("difference of the means: ", s["mu_diff", "median"],
      " [", s["mu_diff", "HDIlo"],
      ", ", s["mu_diff", "HDIup"] , "]\n",sep = "")
  cat("\n")
  cat("The difference of the means is greater than", s["mu_diff","comp"] ,
      "by a probability of", s["mu_diff","%>comp"], "\n")
  cat("and less than", s["mu_diff", "comp"] ,
      "by a probability of", s["mu_diff", "%<comp"], "\n")
  cat("\n")
  invisible(NULL)
}

#' @method summary bayes_paired_wilcox_test
#' @export
summary.bayes_paired_wilcox_test <- function(object, ...) {
  s <- round(object$stats, 3)

  cat("  Data\n")
  cat(object$x_name, ", n = ", length(object$x), "\n", sep = "")
  cat(object$y_name, ", n = ", length(object$y), "\n", sep = "")
  cat("\n")

  #print_bayes_two_sample_t_test_params(object) #Replace by shorter:
  cat("  Model parameters and generated quantities\n")
  cat("mu_diff: the difference in means of ", object$x_name, "and",
      object$y_name, "\n")
  cat("\n")

  cat("  Measures\n" )
  print(s[, c("mean", "sd", "HDIlo", "HDIup", "%<comp", "%>comp")])
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

#adapted from plot.bayes_binom_test
#' @export
plot.bayes_paired_wilcox_test <- function(x, ...) {
  old_par <- par( mar = c(3.5,3.5,2.5,0.5),
                  mgp = c(2.25,0.7,0),
                  mfcol = c(1,1))
  stats <- x$stats
  mcmc_samples <- x$mcmc_samples
  samples_mat <- as.matrix(mcmc_samples)
  mu = samples_mat[,"mu"]
  sample_mat <- as.matrix(x$mcmc_samples)
  xlim = range(c(mu, x$comp))
  plotPost(sample_mat[, "mu_diff"], cred_mass = x$cred_mass,
           comp_val = x$comp, xlim = xlim, cex = 1, cex.lab = 1.5,
           main = "Difference in means",
           xlab = expression(mu_diff), show_median = TRUE)
  par(old_par)
  invisible(NULL)
}

#Adapted from model.code.bayes_binom_test
#output code not yet running!
#' Prints code that replicates the model you just ran.
#'
#' This is good if you better want to understand how the model is
#' implemented or if you want to run a modified version of the code.
#'
#'
#' @param fit The output from a Bayesian First Aid model.
#' @export

# Note: function not recognized as S3, needs to be called by
# model.code.bayes_wilcox_test()
model.code.bayes_paired_wilcox_test <- function(fit) {
  cat("### Model code for the Bayesian First Aid alternative to the binomial",
      "test ###\n\n")
  cat("require(rjags)\n\n")

  cat("# Setting up the data\n")
  cat("x <-", fit$x_name, "\n")
  cat("y <-", fit$y_name, "\n")
  cat("\n")
  pretty_print_function_body(paired_samples_wilcox_model_code)
  invisible(NULL)
}

# Not to be run, just to be printed - adapted from two_sample_t_test_model_code
paired_samples_wilcox_model_code <- function(x, y) {
  d <- NULL
  # The model string written in the JAGS language
  BayesianFirstAid::replace_this_with_model_string

  # Setting parameters for the priors that in practice will result
  # in flat priors on the mu's and sigma's.
  mean_mu = mean( c(x, y), trim = 0.2)
  precision_mu = 1 / (mad( c(x, y) )^2 * 1000000)
  sigma_low = mad( c(x, y) ) / 1000
  sigma_high = mad( c(x, y) ) * 1000

  # Initializing parameters to sensible starting values helps the convergence
  # of the MCMC sampling. Here using robust estimates of the mean (trimmed)
  # and standard deviation (MAD).
  inits_list <- list(mu = mean(d, trim = 0.2), sigma = mad0(d))

  data_list <- list(
    d = d,
    mean_mu = mean(d, trim = 0.2),
    precision_mu = 1 / (mad(d)^2 * 1000000),
    sigma_low = mad(d) / 1000,
    sigma_high = mad(d) * 1000
  )

  # The parameters to monitor.
  params <- c("mu", "sigma")

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
#' @param x The output from a Bayesian First Aid model.
#' @export
# Note: function not recognized as S3, needs to be called by
# model.code.bayes_wilcox_test()

diagnostics.bayes_paired_wilcox_test <- function(x) {
  print_mcmc_info(x$mcmc_samples)
  cat("\n")
  print_diagnostics_measures(round(x$stats, 3))
  cat("\n")
  cat("  Model parameters and generated quantities\n")
  cat("mu_diff: the difference in means of ", x$x_name, "and",
      x$y_name, "\n")
  cat("sigma: the difference in scale of", x$x_name, "and",
      x$y_name, "\n")

  cat("\n")

  old_par <- par( mar = c(3.5,2.5,2.5,0.5) , mgp = c(2.25,0.7,0) )
  plot(x$mcmc_samples)
  par(old_par)
  invisible(NULL)
}

### End wilcox.test S3 functions ###
