#' @export
#'
print.sig_eta_test_bin_wbkg <- function(obj){
  eta_hat <- obj$eta_hat
  se <- obj$SE
  test.stat <- obj$test.stat
  p.val <- obj$p.val
  CI <- obj$CI
  k <- obj$n_bins
  alternative <- obj$alternative
  conf.level <- obj$conf.level
  null.value <- obj$null.value
  beta_hat <- obj$beta.estimate
  M <- obj$samp.sizes[1]; N <- obj$samp.sizes[2]

  cat("Test for signal proportion in the presence of a background-only sample\n\n")

  cat(sprintf("Size of the physics data: %d \nSize of the background-only sample: %d \nNumber of bins: %d",
              N, M, k))
  cat(sprintf("\nstatistic = %.4f, p-value = %.4f\n",
              test.stat, p.val))
  if(alternative == 'two.sided'){
    alt_statement <- "alternative hypothesis: true signal proprotion is not equal to "
  }else if(alternative == 'greater'){
    alt_statement <- "alternative hypothesis: true signal proprotion is greater than "
  }else{
    alt_statement <- "alternative hypothesis: true signal proprotion is less than "
  }

  cat(alt_statement, null.value, "\n", sep = "")

  cat(sprintf("%.0f percent confidence interval of eta:\n",
              conf.level * 100))

  cat(sprintf(" %.4f  %.4f\n",
              CI[1], CI[2]))

  cat("Sample estimates:\n")
  cat(sprintf("Estimate of eta: %.4f, and std.err: %.4f\n", eta_hat, se))
  if(!is.null(beta_hat)){
    cat("Estimate of beta \n")
    cat(sprintf("%.4f\n", beta_hat))
  }
}

#' @export
#'
print.sig_eta_test_unb_wbkg <- function(obj){
  eta_hat <- obj$eta_hat
  se <- obj$SE
  test.stat <- obj$test.stat
  p.val <- obj$p.val
  CI <- obj$CI
  alternative <- obj$alternative
  conf.level <- obj$conf.level
  null.value <- obj$null.value
  beta_hat <- obj$beta.estimate
  m <- obj$samp.sizes[1]; n <- obj$samp.sizes[2]

  cat("Test for signal proportion in the presence of a background-only sample\n\n")

  cat(sprintf("Size of the physics data: %d \nSize of the background-only sample: %d", n, m))

  cat(sprintf("\nstatistic = %.4f, p-value = %.4f\n",
              test.stat, p.val))
  if(alternative == 'two.sided'){
    alt_statement <- "alternative hypothesis: true signal proprotion is not equal to "
  }else if(alternative == 'greater'){
    alt_statement <- "alternative hypothesis: true signal proprotion is greater than "
  }else{
    alt_statement <- "alternative hypothesis: true signal proprotion is less than "
  }

  cat(alt_statement, null.value, "\n", sep = "")

  cat(sprintf("%.0f percent confidence interval of eta:\n",
              conf.level * 100))

  cat(sprintf(" %.4f  %.4f\n",
              CI[1], CI[2]))

  cat("Sample estimates:\n")
  cat(sprintf("Estimate of eta: %.4f, and std.err: %.4f\n", eta_hat, se))
  if(!is.null(beta_hat)){
    cat("Estimate of beta \n")
    cat(sprintf("%.4f\n", beta_hat))
  }
}

#' @export
#'
print.sig_eta_test_bin_wobkg <- function(obj){

  theta0_hat <- obj$theta0_hat
  se <- obj$SE
  test.stat <- obj$test.stat
  p.val <- obj$p.val
  CI <- obj$CI
  k <- obj$bin_counts
  lambda<- obj$lambda
  alternative <- obj$alternative
  conf.level <- obj$conf.level
  null.value <- obj$null.value
  alpha_hat <- obj$alpha.estimate
  N <- obj$samp.size

  cat("Test for signal proportion in the absence of a background-only sample\n\n")

  cat(sprintf("Size of the physics data: %d \n", N))
  cat(sprintf("Number of bins: %d \n", k))
  cat(sprintf("Parameter for the dominating component: %f \n", lambda))
  cat(sprintf("\nstatistic = %.4f, p-value = %.4f\n",
              test.stat, p.val))
  if(alternative == 'two.sided'){
    alt_statement <- "alternative hypothesis: theta0 is not equal to "
  }else if(alternative == 'greater'){
    alt_statement <- "alternative hypothesis: theta0 is greater than "
  }else{
    alt_statement <- "alternative hypothesis: theta0 is less than "
  }

  cat(alt_statement, null.value, "\n", sep = "")

  cat(sprintf("%.0f percent confidence interval of theta0:\n",
              conf.level * 100))

  cat(sprintf(" %.4f  %.4f\n",
              CI[1], CI[2]))

  cat("Sample estimates:\n")
  cat(sprintf("Estimate of theta0 : %.4f, and std.err: %.4f\n", theta0_hat, se))
  cat("Estimate of alpha \n")
  cat(sprintf("%.4f\n", alpha_hat))
}

#' @export
#'
print.sig_eta_test_unb_wobkg <- function(obj){

  theta0_hat <- obj$theta0_hat
  se <- obj$SE
  test.stat <- obj$test.stat
  p.val <- obj$p.val
  CI <- obj$CI
  lambda<- obj$lambda
  alternative <- obj$alternative
  conf.level <- obj$conf.level
  null.value <- obj$null.value
  alpha_hat <- obj$alpha.estimate
  n <- obj$samp.size

  cat("Test for signal proportion in the absence of a background-only sample\n\n")

  cat(sprintf("Size of the physics data: %d \n", n))
  cat(sprintf("Parameter for the dominating component: %f \n", lambda))
  cat(sprintf("\nstatistic = %.4f, p-value = %.4f\n",
              test.stat, p.val))
  if(alternative == 'two.sided'){
    alt_statement <- "alternative hypothesis: theta0 is not equal to "
  }else if(alternative == 'greater'){
    alt_statement <- "alternative hypothesis: theta0 is greater than "
  }else{
    alt_statement <- "alternative hypothesis: theta0 is less than "
  }

  cat(alt_statement, null.value, "\n", sep = "")

  cat(sprintf("%.0f percent confidence interval of theta0:\n",
              conf.level * 100))

  cat(sprintf(" %.4f  %.4f\n",
              CI[1], CI[2]))

  cat("Sample estimates:\n")
  cat(sprintf("Estimate of theta0 : %.4f, and std.err: %.4f\n", theta0_hat, se))
  cat("Estimate of alpha \n")
  cat(sprintf("%.4f\n", alpha_hat))
}
