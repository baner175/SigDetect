sig_detect_unb_par_Yden_wobkg <- function(phys.samp,
                                          search_region,
                                          fs, bench.bkg,
                                          d_log_bkg, d2_log_bkg,
                                          eta0,
                                          lambda,
                                          mu1, mu2, sig0,
                                          alpha_lower,
                                          alpha_upper,
                                          alpha_start,
                                          conf.level,
                                          alternative,
                                          CI_type){

  n <- length(phys.samp)
  l <- search_region[1]; u <- search_region[2]
  p <- length(alpha_start)
  neg_loglik <- function(a){
    qi <- sapply(phys.samp, bench.bkg, alpha = a)
    return(-sum(log(qi)))
  }
  alpha_hat <- nlminb(start = alpha_start,
                      objective = neg_loglik,
                      lower = alpha_lower,
                      upper = alpha_upper)$par
  q <- function(x) bench.bkg(x, alpha = alpha_hat)
  g <- function(x){
    phi1 <- truncdist::dtrunc(x, spec = 'norm', a = l, b = u,
                   mean = mu1, sd = sig0)
    phi2 <- truncdist::dtrunc(x, spec = 'norm', a = l, b = u,
                   mean = mu2, sd = sig0)
    return(
      (1-2*lambda)*q(x) + lambda*(phi1+phi2)
    )
  }
  S <- function(x){
    return(fs(x)/g(x) -1)
  }
  normS <- integrate(function(x){
    (S(x)^2)*g(x)
  }, l,  u)$value |> sqrt()
  S0 <- function(x) {
    S(x)/(normS^2)
  }

  if(length(formalArgs(d_log_bkg)) == 2){
    d_log_q_hat <- function(x) {
      val <- d_log_bkg(x, alpha = alpha_hat)
      dim(val) <- c(p,1)
      return(val)
    }
  }else{
    d_log_q_hat <- function(x) {
      val <- d_log_bkg(x)
      dim(val) <- c(p,1)
      return(val)
    }
  }

  if(class(d2_log_bkg)[1] == 'function'){
    if(length(formalArgs(d2_log_bkg)) == 2){
      d2_log_q_hat <- function(x){
        val <- d2_log_bkg(x, alpha = alpha_hat)
        dim(val) <- c(p,p)
        return(val)
      }
    }else{
      d2_log_q_hat <- d2_log_bkg(alpha_hat)
      dim(d2_log_q_hat) <- c(p,p)
    }
  }else{
    d2_log_q_hat <- d2_log_bkg
  }

  d_normS_sq <- sapply(1:p, function(i){
    -(1-2*lambda)*integrate(function(x) {
      sapply(x, function(xi) {
        fs_val <- fs(xi)
        g_val  <- g(xi)
        q_val <- q(xi)
        d_log_q_val <- d_log_q_hat(xi)[i]
        return(((fs_val/g_val)^2)*q_val*d_log_q_val)
      })
    }, l, u)$value
  })
  dim(d_normS_sq) <- c(p,1)

  d_S0 <- function(x){
    num <- (normS^2)*(S(x)+1)*(q(x)/g(x))*(1-2*lambda)*d_log_q_hat(x) + S(x)*d_normS_sq
    return(-num/(normS^4))
  }
  S0_phys <- sapply(phys.samp, S0)

  theta0_hat <- mean(S0_phys)
  test_num <- sqrt(n)*(theta0_hat - eta0)

  var_S0_F <- mean(sapply(phys.samp, S0)^2) - theta0_hat^2
  d_S0_phys <- sapply(phys.samp, d_S0)
  dim(d_S0_phys) <- c(p,n)
  d_theta0 <- apply(d_S0_phys, MARGIN = 1, mean)
  dim(d_theta0) <- c(p,1)

  if(class(d2_log_q_hat)[1] == 'function'){
    d2_log_q_list <- lapply(phys.samp, d2_log_q_hat)
    J_hat <- -Reduce("+", d2_log_q_list)/n
  }else{
    J_hat <- -d2_log_q_hat
  }
  J_hat_inv <- solve(J_hat)
  V_hat_list <- lapply(phys.samp, function(x){
    val <- d_log_q_hat(x)
    dim(val) <- c(p,1)
    return(val %*% t(val))
  })
  V_hat <- Reduce("+", V_hat_list)/n

  Cov_integrand <- function(x) {
    val <- d_log_q_hat(x)*S0(x)
    dim(val) <- c(p,1)
    return(val)
  }
  Cov_integrand_phys <- sapply(phys.samp, Cov_integrand)
  dim(Cov_integrand_phys) <- c(p,n)
  C_hat <- apply(Cov_integrand_phys, 1, mean)
  test_denom <- sqrt(
    var_S0_F +
      as.numeric(t(d_theta0)%*%J_hat_inv%*%V_hat%*%J_hat_inv%*%d_theta0) +
      2*as.numeric(t(d_theta0)%*%J_hat_inv%*%C_hat)
  )

  std_err <- test_denom/sqrt(n)
  test.stat <- test_num/test_denom
  alpha_val <- 1-conf.level
  if(alternative == 'greater'){
    p.val <- pnorm(test.stat, lower.tail = FALSE)
  }else if(alternative == 'lower'){
    p.val <- pnorm(test.stat)
  }else{
    p.val <- 2*pnorm(-abs(test.stat))
  }
  if(CI_type == 'greater'){
    CI <- c(theta0_hat - qnorm(alpha_val)*std_err, Inf)
  }else if(CI_type == 'lower'){
    CI <- c(-Inf, theta0_hat + qnorm(alpha_val, lower.tail = FALSE)*std_err)
  }else{
    CI <- theta0_hat + c(-1,1)*qnorm(alpha_val/2, lower.tail = FALSE)*std_err
  }
  return(list(theta0_hat = theta0_hat, SE = std_err, test.stat = test.stat, p.val = p.val, CI = CI,
              samp.size = n, alternative = alternative, lambda = lambda,
              conf.level = conf.level, null.value = eta0, alpha.estimate = alpha_hat))
}

sig_detect_unb_par_Nden_wobkg <- function(phys.samp,
                                          search_region,
                                          fs, bench.bkg,
                                          d_log_bkg, d2_log_bkg,
                                          eta0,
                                          lambda,
                                          mu1, mu2, sig0,
                                          alpha_lower,
                                          alpha_upper,
                                          alpha_start,
                                          conf.level,
                                          alternative,
                                          CI_type){
  n <- length(phys.samp)
  l <- search_region[1]; u <- search_region[2]
  p <- length(alpha_start)

  neg_loglik <- function(a){
    h_mass <- integrate(bench.bkg, l, u, alpha = a)$value
    qi <- sapply(phys.samp, bench.bkg, alpha = a)/h_mass
    return(-sum(log(qi)))
  }
  alpha_hat <- nlminb(start = alpha_start,
                      objective = neg_loglik,
                      lower = alpha_lower,
                      upper = alpha_upper)$par
  h_mass <- integrate(bench.bkg, l, u, alpha = alpha_hat)$value
  h <- function(x) bench.bkg(x, alpha = alpha_hat)
  q <- function(x) bench.bkg(x, alpha = alpha_hat)/h_mass
  g <- function(x){
    phi1 <- truncdist::dtrunc(x, spec = 'norm', a = l, b = u,
                   mean = mu1, sd = sig0)
    phi2 <- truncdist::dtrunc(x, spec = 'norm', a = l, b = u,
                   mean = mu2, sd = sig0)
    return((1-2*lambda)*q(x) + lambda*(phi1+phi2))
  }

  if(length(formalArgs(d_log_bkg)) == 2){
    d_log_h_hat <- function(x) {
      val <- d_log_bkg(x, alpha = alpha_hat)
      dim(val) <- c(p,1)
      return(val)
    }
  }else{
    d_log_h_hat <- function(x){
      val <- d_log_bkg(x)
      dim(val) <- c(p,1)
      return(val)
    }
  }

  S <- function(x){
    return(fs(x)/g(x) -1)
  }
  normS <- integrate(function(x){
    (S(x)^2)*g(x)
  }, l,  u)$value |> sqrt()
  S0 <- function(x) {
    S(x)/(normS^2)
  }
  E_q_d_log_h <- sapply(1:p, function(i){
    integrate(function(x){
      sapply(x, function(xi){
        q(xi)*d_log_h_hat(xi)[i,1]
      })
    }, l, u)$value
  })
  dim(E_q_d_log_h) <- c(p,1)
  d_log_q_hat <- function(x) {
    val <- d_log_h_hat(x) - E_q_d_log_h
    dim(val) <- c(p,1)
    return(val)
  }


  if(class(d2_log_bkg)[1] == 'function'){
    if(length(formalArgs(d2_log_bkg)) == 2){
      d2_log_h_hat <- function(x){
        val <- d2_log_bkg(x, alpha = alpha_hat)
        dim(val) <- c(p,p)
        return(val)
      }
    }else{
      d2_log_h_hat <- d2_log_bkg(alpha_hat)
      dim(d2_log_h_hat) <- c(p,p)
    }
  }else{
    d2_log_h_hat <- d2_log_bkg
  }

  d2_mat2_fun <- function(x){
    d_log_h_val <- d_log_h_hat(x)
    return(d_log_h_val %*% t(d_log_h_val) * q(x))
  }
  d2_mat3 <- E_q_d_log_h %*% t(E_q_d_log_h)
  if(class(d2_log_h_hat)[1] == 'function'){
    d2_mat1 <- matrix(0,p,p)
    d2_mat2 <- matrix(0, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        d2_mat1[i, j] <- integrate(function(x) {
          sapply(x, function(xi) q(xi)*d2_log_h_hat(xi)[i, j])
        }, l, u)$value
        d2_mat2[i, j] <- integrate(function(x) {
          sapply(x, function(xi) d2_mat2_fun(xi)[i, j])
        }, l, u)$value
      }
    }
    d2_log_q_hat <- function(x){
      d2_log_h_hat(x) - d2_mat1 - d2_mat2 + d2_mat3
    }
    d2_log_q_list <- lapply(phys.samp, d2_log_q_hat)
    J_hat <- -Reduce("+", d2_log_q_list)/n
  }else{
    d2_mat1 <- d2_log_h_hat
    d2_mat2 <- matrix(0, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        d2_mat2[i, j] <- integrate(function(x) {
          sapply(x, function(xi) d2_mat2_fun(xi)[i, j])
        }, l, u)$value
      }
    }
    d2_log_q_hat <- d2_log_h_hat - d2_mat1 - d2_mat2 + d2_mat3
    J_hat <- -d2_log_q_hat
  }


  d_normS_sq <- sapply(1:p, function(i){
    -(1-2*lambda)*integrate(function(x) {
      sapply(x, function(xi) {
        fs_val <- fs(xi)
        g_val  <- g(xi)
        q_val <- q(xi)
        d_log_q_val <- d_log_q_hat(xi)[i]
        return(((fs_val/g_val)^2)*q_val*d_log_q_val)
      })
    }, l, u)$value
  })
  dim(d_normS_sq) <- c(p,1)

  d_S0 <- function(x){
    num <- (normS^2)*(S(x)+1)*(q(x)/g(x))*d_log_q_hat(x)*(1-2*lambda) + S(x)*d_normS_sq
    return(-num/(normS^4))
  }

  S0_phys <- sapply(phys.samp, S0)

  theta0_hat <- mean(S0_phys)
  test_num <- sqrt(n)*(theta0_hat - eta0)

  var_S0_F <- mean(sapply(phys.samp, S0)^2) - theta0_hat^2
  d_S0_phys <- sapply(phys.samp, d_S0)
  dim(d_S0_phys) <- c(p,n)
  d_theta0 <- apply(d_S0_phys, MARGIN = 1, mean)
  dim(d_theta0) <- c(p,1)


  J_hat_inv <- solve(J_hat)
  V_hat_list <- lapply(phys.samp, function(x){
    val <- d_log_q_hat(x)
    dim(val) <- c(p,1)
    return(val %*% t(val))
  })
  V_hat <- Reduce("+", V_hat_list)/n

  Cov_integrand <- function(x) {
    val <- d_log_q_hat(x)*S0(x)
    dim(val) <- c(p,1)
    return(val)
  }
  Cov_integrand_phys <- sapply(phys.samp, Cov_integrand)
  dim(Cov_integrand_phys) <- c(p,n)
  C_hat <- apply(Cov_integrand_phys, 1, mean)
  test_denom <- sqrt(
    var_S0_F +
      as.numeric(t(d_theta0)%*%J_hat_inv%*%V_hat%*%J_hat_inv%*%d_theta0) +
      2*as.numeric(t(d_theta0)%*%J_hat_inv%*%C_hat)
  )
  std_err <- test_denom/sqrt(n)
  test.stat <- test_num/test_denom
  alpha_val <- 1-conf.level

  if(alternative == 'greater'){
    p.val <- pnorm(test.stat, lower.tail = FALSE)
  }else if(alternative == 'lower'){
    p.val <- pnorm(test.stat)
  }else{
    p.val <- 2*pnorm(-abs(test.stat))
  }
  if(CI_type == 'greater'){
    CI <- c(theta0_hat - qnorm(alpha_val)*std_err, Inf)
  }else if(CI_type == 'lower'){
    CI <- c(-Inf, theta0_hat + qnorm(alpha_val, lower.tail = FALSE)*std_err)
  }else{
    CI <- theta0_hat + c(-1,1)*qnorm(alpha_val/2, lower.tail = FALSE)*std_err
  }
  return(list(theta0_hat = theta0_hat, SE = std_err, test.stat = test.stat, p.val = p.val, CI = CI,
              samp.size = n, alternative = alternative, lambda = lambda,
              conf.level = conf.level, null.value = eta0, alpha.estimate = alpha_hat))
}


#' Performs a conservative test on the signal proportion using unbinned data data in the absence of a background-only dataset
#'
#' @param phys.samp a numeric vector containing the physics sample.
#' @param search_region a numeric vector of two numbers with the first and the second number being the lower and upper limits of the search region respectively.
#' @param fs the function for the density of the signal distribution.
#' @param bench.bkg a function specifying the benchmark background model with arguments `x` and `alpha`. The function is not required to be a density.
#' @param d_log_bkg a function that computes the score function of the `bench.bkg` that may or may not depend on the parameter `alpha`.
#' @param d2_log_bkg a function that computes the Hessian matrix of the `bench.bkg` with respect to `alpha` or a numeric value.
#' @param eta0 a numeric value for the null value of the signal proportion.
#' @param is.density a logical indicating whether `bench.bkg` is a density.
#' @param lambda a numeric value determining the size of the diffused bump added to `bench.bkg`. The value must be between 0 and 0.5.
#' @param mu1 a numeric value for the location of the first Gaussian component in the diffused bump.
#' @param mu2 a numeric value for the location of the second Gaussian component in the diffused bump.
#' @param sig0 a numeric value for the common scale parameter of the Gaussian components in the diffused bump.
#' @param alpha_lower a numeric vector or a numeric value for the lower limit of `alpha`, the parameter in `bench.bkg`.
#' @param alpha_upper a numeric vector or a numeric value for the upper limit of `alpha`, the parameter in `bench.bkg`.
#' @param alpha_start a numeric vector or a numeric value for the initial value of `alpha` to obtains its MLE.
#' @param conf.level confidence level of the interval.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param CI_type a character string specifying the type of confidence interval to compute, must be one of "two.sided" (default), "greater" or "less".
#' @return a list with class "sig_eta_test_unb_wobkg" containing the following components:
#' \describe{
#'   \item{theta0_hat}{estimated theta0, a conservative estimate to the true signal proportion.}
#'   \item{SE}{standard error for `theta0_hat`}
#'   \item{test.stat}{observed test statistic.}
#'   \item{p.val}{the p-value for the test.}
#'   \item{CI}{a confidence interval for theta0.}
#'   \item{samp.size}{size of the physics sample.}
#'   \item{alternative}{a character string describing the alternative hypothesis.}
#'   \item{lambda}{the size of the diffused bump used for testing.}
#'   \item{conf.level}{confidence level for the confidence interval.}
#'   \item{null.value}{the specified hypothesized value of the signal proportion.}
#'   \item{alpha.estimate}{MLE of the parameter in `bench.bkg`, `alpha`.}
#' }
#'
#' @export

sig_detect_unbinned_wobkg <- function(phys.samp,
                                      search_region,
                                      fs, bench.bkg,
                                      d_log_bkg = NULL, d2_log_bkg = NULL,
                                      eta0 = 0,
                                      is.density = TRUE,
                                      lambda = 0,
                                      mu1, mu2, sig0,
                                      alpha_lower = NULL,
                                      alpha_upper = NULL,
                                      alpha_start = NULL,
                                      conf.level = 0.95,
                                      alternative = c("greater", "two.sided", "less"),
                                      CI_type = c("two.sided", "greater", "less")){

  arg_names <- formalArgs(bench.bkg)
  n <- length(phys.samp)
  l <- search_region[1]; u <- search_region[2]
  lambda <- lambda
  p <- length(alpha_start)
  alternative <- match.arg(alternative)
  CI_type <- match.arg(CI_type)
  if(length(formalArgs(fs))>1) stop('fs must only have one argument')
  if((lambda < 0) | (lambda > 0.5)) stop('lambda must be between 0 and 0.5')
  if(length(arg_names) != 2) stop('bench.bkg must have two arguments, x and alpha')
  if((length(arg_names) == 2) & ((arg_names[1] != 'x') | (arg_names[2] != 'alpha'))) stop('Arguments of bench.bkg must be x and alpha')
  if((length(arg_names) == 2) & is.null(d_log_bkg)) stop('Must provide the gradient vector d_log_bkg')
  if((length(arg_names) == 2) & is.null(d2_log_bkg)) stop('Must provide the Hessian matrix d2_log_bkg')
  if((length(arg_names) == 2) & is.null(alpha_lower)) stop('Must provide the lower limits for alpha')
  if((length(arg_names) == 2) & is.null(alpha_upper)) stop('Must provide the upper limits for alpha')
  if((length(arg_names) == 2) & is.null(alpha_start)) stop('Must provide an initial value for alpha to obtain MLE')

  arg_names_d_log_bkg <- formalArgs(d_log_bkg)
  if(length(arg_names_d_log_bkg) > 2) stop('d_log_bkg must have at most two arguments')
  if((length(arg_names_d_log_bkg) == 2) & (arg_names_d_log_bkg[2] != 'alpha')) stop('The second argument of d_log_bkg must be alpha')

  if(class(d2_log_bkg)[1] == 'function'){
    arg_names_d2_log_bkg <- formalArgs(d2_log_bkg)
    if(length(arg_names_d2_log_bkg) > 2) stop('d2_log_bkg must have at most two arguments, x and alpha')
    if((length(arg_names_d2_log_bkg) == 2) & (arg_names_d2_log_bkg[2] != 'alpha')) stop('The second argument of d2_log_bkg must be alpha')
  }else if(d2_log_bkg != 0){
    stop("The argument d2_log_bkg must be a function or the numeric value 0")
  }else{
    d2_log_bkg <- matrix(0, p, p)
  }

  if(is.density){
    res <- sig_detect_unb_par_Yden_wobkg(phys.samp = phys.samp,
                                         search_region = c(l,u),
                                         fs = fs, bench.bkg = bench.bkg,
                                         d_log_bkg = d_log_bkg, d2_log_bkg = d2_log_bkg,
                                         eta0 = eta0,
                                         lambda = lambda,
                                         mu1 = mu1, mu2 = mu2, sig0 = sig0,
                                         alpha_lower = alpha_lower,
                                         alpha_upper = alpha_upper,
                                         alpha_start = alpha_start,
                                         conf.level = conf.level,
                                         alternative = alternative,
                                         CI_type = CI_type)
    class(res) <- 'sig_eta_test_unb_wobkg'
    return(res)

  }else{
    res <- sig_detect_unb_par_Nden_wobkg(phys.samp = phys.samp,
                                         search_region = c(l,u),
                                         fs = fs, bench.bkg = bench.bkg,
                                         d_log_bkg = d_log_bkg, d2_log_bkg = d2_log_bkg,
                                         eta0 = eta0,
                                         lambda = lambda,
                                         mu1 = mu1, mu2 = mu2, sig0 = sig0,
                                         alpha_lower = alpha_lower,
                                         alpha_upper = alpha_upper,
                                         alpha_start = alpha_start,
                                         conf.level = conf.level,
                                         alternative = alternative,
                                         CI_type = CI_type)
    class(res) <- 'sig_eta_test_unb_wobkg'
    return(res)
  }
}
