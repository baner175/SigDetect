sig_detect_bin_par_Yden <- function(phys.counts,
                                    bkg.counts,
                                    search_region,
                                    fs, prop.bkg,
                                    d_log_bkg, d2_log_bkg,
                                    eta0,
                                    beta_lower,
                                    beta_upper,
                                    beta_start,
                                    conf.level,
                                    alternative,
                                    CI_type){

  mi <- bkg.counts; ni <- phys.counts
  M <- sum(mi); N <- sum(ni)
  k <- length(bkg.counts)
  l <- search_region[1]; u <- search_region[2]
  bin_ends <- seq(l, u, length.out = k+1)
  xi <- (bin_ends[1:k] + bin_ends[2:(k+1)])/2
  p <- length(beta_start)
  neg_loglik <- function(b){
    gi <- sapply(1:k, function(i){
      integrate(prop.bkg, lower = bin_ends[i], upper = bin_ends[i+1], beta = b)$value
    })
    return(-sum(mi*log(gi)))
  }
  beta_hat <- nlminb(start = beta_start,
                     objective = neg_loglik,
                     lower = beta_lower,
                     upper = beta_upper)$par
  g <- function(x) prop.bkg(x, beta = beta_hat)
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
    d_log_g_hat <- function(x) {
      val <- d_log_bkg(x, beta = beta_hat)
      dim(val) <- c(p,1)
      return(val)
    }
  }else{
    d_log_g_hat <- function(x) {
      val <- d_log_bkg(x)
      dim(val) <- c(p,1)
      return(val)
    }
  }

  if(class(d2_log_bkg)[1] == 'function'){
    if(length(formalArgs(d2_log_bkg)) == 2){
      d2_log_g_hat <- function(x){
        val <- d2_log_bkg(x, beta = beta_hat)
        dim(val) <- c(p,p)
        return(val)
      }
    }else{
      d2_log_g_hat <- d2_log_bkg(beta_hat)
      dim(d2_log_g_hat) <- c(p,p)
    }
  }else{
    d2_log_g_hat <- d2_log_bkg
  }

  d_normS_sq <- sapply(1:p, function(i){
    -integrate(function(x) {
      sapply(x, function(xi) {
        fs_val <- fs(xi)
        g_val <- g(xi)
        d_log_g_val <- d_log_g_hat(xi)[i]
        ((fs_val^2) / g_val) * d_log_g_val
      })
    }, l, u)$value
  })
  dim(d_normS_sq) <- c(p,1)

  d_S0 <- function(x){
    num <- (normS^2)*(S(x)+1)*d_log_g_hat(x) + S(x)*d_normS_sq
    return(-num/(normS^4))
  }

  S0_xi <- sapply(xi, S0)

  theta0_hat <- sum(S0_xi*ni)/N; delta0_hat <- sum(S0_xi*mi)/M
  eta_hat <- (theta0_hat - delta0_hat)/(1-delta0_hat)
  test_num <- sqrt(M*N)*(eta_hat - eta0)

  sig_theta0_sq <- sum((S0_xi^2)*ni)/N - theta0_hat^2
  sig_delta0_sq <- sum((S0_xi^2)*mi)/M - delta0_hat^2
  dT_1 <- 1/(1-delta0_hat); dT_2 <- (theta0_hat - 1)/((1-delta0_hat)^2)
  d_S0_xi <- sapply(xi, d_S0)
  dim(d_S0_xi) <- c(p,k)
  d_theta0 <- (d_S0_xi %*% ni)/N
  d_delta0 <- (d_S0_xi %*% mi)/M
  Gamma_hat <- dT_1*d_theta0 + dT_2*d_delta0

  if(class(d2_log_g_hat)[1] == 'function'){
    d2_log_g_list <- lapply(xi, d2_log_g_hat)
    temp_mat <- matrix(0,p,p)
    for(i in 1:k){
      temp_mat <- temp_mat + d2_log_g_list[[i]]*mi[i]
    }
    J_hat <- -temp_mat/M
  }else{
    J_hat <- -d2_log_g_hat
  }
  J_hat_inv <- solve(J_hat)
  V_hat_list <- lapply(xi, function(x){
    val <- d_log_g_hat(x)
    dim(val) <- c(p,1)
    return(val %*% t(val))
  })
  temp_mat <- matrix(0,p,p)
  for(i in 1:k){
    temp_mat <- temp_mat + V_hat_list[[i]]*mi[i]
  }
  V_hat <- temp_mat/M

  Cov_integrand <- function(x) {
    val <- d_log_g_hat(x)*S0(x)
    dim(val) <- c(p,1)
    return(val)
  }
  Cov_integrand_xi <- sapply(xi, Cov_integrand)
  dim(Cov_integrand_xi) <- c(p,k)
  C_hat <- (Cov_integrand_xi %*% mi)/M
  test_denom <- (M*(dT_1^2)*sig_theta0_sq + N*(
    (dT_2^2)*sig_delta0_sq + as.numeric(t(Gamma_hat)%*%J_hat_inv%*%V_hat%*%J_hat_inv%*%Gamma_hat) +
      2*dT_2*as.numeric(t(Gamma_hat)%*%J_hat_inv%*%C_hat)
  )) |> sqrt()
  sig_eta <- test_denom/sqrt(M+N)
  std_err <- sig_eta*sqrt(M+N)/sqrt(M*N)
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
    CI <- c(eta_hat - qnorm(alpha_val)*std_err, Inf)
  }else if(CI_type == 'lower'){
    CI <- c(-Inf, eta_hat + qnorm(alpha_val, lower.tail = FALSE)*std_err)
  }else{
    CI <- eta_hat + c(-1,1)*qnorm(alpha_val/2, lower.tail = FALSE)*std_err
  }
  return(list(eta_hat = eta_hat, SE = std_err, test.stat = test.stat, p.val = p.val, CI = CI,
              samp.sizes = c(M,N), n_bins = k, alternative = alternative,
              conf.level = conf.level, null.value = eta0, beta.estimate = beta_hat))
}

sig_detect_bin_par_Nden <- function(phys.counts,
                                    bkg.counts,
                                    search_region,
                                    fs, prop.bkg,
                                    d_log_bkg, d2_log_bkg,
                                    eta0,
                                    beta_lower,
                                    beta_upper,
                                    beta_start,
                                    conf.level,
                                    alternative,
                                    CI_type){

  mi <- bkg.counts; ni <- phys.counts
  M <- sum(mi); N <- sum(ni)
  k <- length(bkg.counts)
  l <- search_region[1]; u <- search_region[2]
  bin_ends <- seq(l, u, length.out = k+1)
  xi <- (bin_ends[1:k] + bin_ends[2:(k+1)])/2
  p <- length(beta_start)

  neg_loglik <- function(b){
    h_mass <- integrate(prop.bkg, l, u, beta = b)$value
    gi <- sapply(1:k, function(i){
      integrate(prop.bkg, lower = bin_ends[i], upper = bin_ends[i+1], beta = b)$value/h_mass
    })
    return(-sum(mi*log(gi)))
  }
  beta_hat <- nlminb(start = beta_start,
                     objective = neg_loglik,
                     lower = beta_lower,
                     upper = beta_upper)$par
  h_mass <- integrate(prop.bkg, l, u, beta = beta_hat)$value
  h <- function(x) prop.bkg(x, beta = beta_hat)
  g <- function(x) prop.bkg(x, beta = beta_hat)/h_mass

  if(length(formalArgs(d_log_bkg)) == 2){
    d_log_h_hat <- function(x) {
      val <- d_log_bkg(x, beta = beta_hat)
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
  E_g_d_log_h <- sapply(1:p, function(i){
    integrate(function(x){
      sapply(x, function(xi){
        g(xi)*d_log_h_hat(xi)[i,1]
      })
    }, l, u)$value
  })
  dim(E_g_d_log_h) <- c(p,1)
  d_log_g_hat <- function(x) {
    val <- d_log_h_hat(x) - E_g_d_log_h
    dim(val) <- c(p,1)
    return(val)
  }

  if(class(d2_log_bkg)[1] == 'function'){
    if(length(formalArgs(d2_log_bkg)) == 2){
      d2_log_h_hat <- function(x){
        val <- d2_log_bkg(x, beta = beta_hat)
        dim(val) <- c(p,p)
        return(val)
      }
    }else{
      d2_log_h_hat <- d2_log_bkg(beta_hat)
      dim(d2_log_h_hat) <- c(p,p)
    }
  }else{
    d2_log_h_hat <- d2_log_bkg
  }

  d2_mat2_fun <- function(x){
    d_log_h_val <- d_log_h_hat(x)
    return(d_log_h_val %*% t(d_log_h_val) * g(x))
  }
  d2_mat3 <- E_g_d_log_h %*% t(E_g_d_log_h)

  if(class(d2_log_h_hat)[1] == 'function'){
    d2_mat1 <- matrix(0,p,p)
    d2_mat2 <- matrix(0, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        d2_mat1[i, j] <- integrate(function(x) {
          sapply(x, function(xi) g(xi)*d2_log_h_hat(xi)[i, j])
        }, l, u)$value
        d2_mat2[i, j] <- integrate(function(x) {
          sapply(x, function(xi) d2_mat2_fun(xi)[i, j])
        }, l, u)$value
      }
    }
    d2_log_g_hat <- function(x){
      d2_log_h_hat(x) - d2_mat1 - d2_mat2 + d2_mat3
    }
    d2_log_g_list <- lapply(xi, d2_log_g_hat)
    temp_mat <- matrix(0,p,p)
    for(i in 1:k){
      temp_mat <- temp_mat + d2_log_g_list[[i]]*mi[i]
    }
    J_hat <- -temp_mat/M
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
    d2_log_g_hat <- - d2_mat2 + d2_mat3
    J_hat <- -d2_log_g_hat
  }


  d_normS_sq <- sapply(1:p, function(i){
    -integrate(function(x) {
      sapply(x, function(xi) {
        fs_val <- fs(xi)
        g_val  <- g(xi)
        d_log_g_val <- d_log_g_hat(xi)[i]
        ((fs_val^2) / g_val) * d_log_g_val
      })
    }, l, u)$value
  })
  dim(d_normS_sq) <- c(p,1)

  d_S0 <- function(x){
    num <- (normS^2)*(S(x)+1)*d_log_g_hat(x) + S(x)*d_normS_sq
    return(-num/(normS^4))
  }

  S0_xi <- sapply(xi, S0)

  theta0_hat <- sum(S0_xi*ni)/N; delta0_hat <- sum(S0_xi*mi)/M
  eta_hat <- (theta0_hat - delta0_hat)/(1-delta0_hat)
  test_num <- sqrt(M*N)*(eta_hat - eta0)

  sig_theta0_sq <- sum((S0_xi^2)*ni)/N - theta0_hat^2
  sig_delta0_sq <- sum((S0_xi^2)*mi)/M - delta0_hat^2
  dT_1 <- 1/(1-delta0_hat); dT_2 <- (theta0_hat - 1)/((1-delta0_hat)^2)
  d_S0_xi <- sapply(xi, d_S0)
  dim(d_S0_xi) <- c(p,k)
  d_theta0 <- (d_S0_xi %*% ni)/N
  d_delta0 <- (d_S0_xi %*% mi)/M
  Gamma_hat <- dT_1*d_theta0 + dT_2*d_delta0
  J_hat_inv <- solve(J_hat)
  V_hat_list <- lapply(xi, function(x){
    val <- d_log_g_hat(x)
    dim(val) <- c(p,1)
    return(val %*% t(val))
  })
  temp_mat <- matrix(0,p,p)
  for(i in 1:k){
    temp_mat <- temp_mat + V_hat_list[[i]]*mi[i]
  }
  V_hat <- temp_mat/M

  Cov_integrand <- function(x) {
    val <- d_log_g_hat(x)*S0(x)
    dim(val) <- c(p,1)
    return(val)
  }
  Cov_integrand_xi <- sapply(xi, Cov_integrand)
  dim(Cov_integrand_xi) <- c(p,k)
  C_hat <- (Cov_integrand_xi %*% mi)/M
  test_denom <- (M*(dT_1^2)*sig_theta0_sq + N*(
    (dT_2^2)*sig_delta0_sq + as.numeric(t(Gamma_hat)%*%J_hat_inv%*%V_hat%*%J_hat_inv%*%Gamma_hat) +
      2*dT_2*as.numeric(t(Gamma_hat)%*%J_hat_inv%*%C_hat)
  )) |> sqrt()
  sig_eta <- test_denom/sqrt(M+N)
  std_err <- sig_eta*sqrt(M+N)/sqrt(M*N)
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
    CI <- c(eta_hat - qnorm(alpha_val)*std_err, Inf)
  }else if(CI_type == 'lower'){
    CI <- c(-Inf, eta_hat + qnorm(alpha_val, lower.tail = FALSE)*std_err)
  }else{
    CI <- eta_hat + c(-1,1)*qnorm(alpha_val/2, lower.tail = FALSE)*std_err
  }
  return(list(eta_hat = eta_hat, SE = std_err, test.stat = test.stat, p.val = p.val, CI = CI,
              samp.sizes = c(M,N), n_bins = k, alternative = alternative,
              conf.level = conf.level, null.value = eta0, beta.estimate = beta_hat))
}

#' Performs a conservative test on the signal proportion using binned data in the presence of a background-only dataset
#'
#' @param phys.counts a numeric vector containing the bin counts for the physics sample.
#' @param bkg.counts a numeric vector containing the bin counts for the background-only sample.
#' @param search_region a numeric vector of two numbers with the first and the second number being the lower and upper limits of the search region respectively.
#' @param fs the function for the density of the signal distribution.
#' @param prop.bkg a function specifying the proposal background model with arguments `x` and `beta`. The function is not required to be a density.
#' @param d_log_bkg a function that computes the score function of the `prop.bkg` that may or may not depend on the parameter `beta`.
#' @param d2_log_bkg a function that computes the Hessian matrix of the `prop.bkg` with respect to `beta` or a numeric value.
#' @param eta0 a numeric value for the null value of the signal proportion.
#' @param is.density a logical indicating whether `prop.bkg` is a density.
#' @param beta_lower a numeric vector or a numeric value for the lower limit of `beta`, the parameter in `bench.bkg`.
#' @param beta_upper a numeric vector or a numeric value for the upper limit of `beta`, the parameter in `bench.bkg`.
#' @param beta_start a numeric vector or a numeric value for the initial value of `beta` to obtains its MLE.
#' @param conf.level confidence level of the interval.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param CI_type a character string specifying the type of confidence interval to compute, must be one of "two.sided" (default), "greater" or "less".
#' @return a list with class "sig_eta_test_bin_wbkg" containing the following components:
#' \describe{
#'   \item{eta_hat}{estimated eta, an estimate to the true signal proportion.}
#'   \item{SE}{standard error for `eta_hat`}
#'   \item{test.stat}{observed test statistic.}
#'   \item{p.val}{the p-value for the test.}
#'   \item{CI}{a confidence interval for eta.}
#'   \item{samp.sizes}{a numeric vector containing the total sizes of the background-only sample and the physics sample.}
#'   \item{n_bins}{number of bins used.}
#'   \item{alternative}{a character string describing the alternative hypothesis.}
#'   \item{conf.level}{confidence level for the confidence interval.}
#'   \item{null.value}{the specified hypothesized value of the signal proportion.}
#'   \item{beta.estimate}{MLE of the parameter in `prop.bkg`, `beta`.}
#' }
#'
#' @export

sig_detect_binned_wbkg <- function(phys.counts,
                                   bkg.counts,
                                   search_region,
                                   fs, prop.bkg,
                                   d_log_bkg = NULL, d2_log_bkg = NULL,
                                   eta0 = 0,
                                   is.density = TRUE,
                                   beta_lower = NULL,
                                   beta_upper = NULL,
                                   beta_start = NULL,
                                   conf.level = 0.95,
                                   alternative = c("greater", "two.sided", "less"),
                                   CI_type = c("two.sided", "greater", "less")){

  mi <- bkg.counts; ni <- phys.counts
  M <- sum(mi); N <- sum(ni)
  k <- length(bkg.counts)
  l <- search_region[1]; u <- search_region[2]
  bin_ends <- seq(l, u, length.out = k+1)
  xi <- (bin_ends[1:k] + bin_ends[2:(k+1)])/2
  p <- length(beta_start)
  alternative <- match.arg(alternative)
  CI_type <- match.arg(CI_type)

  arg_names <- formalArgs(prop.bkg)
  if(length(formalArgs(fs))>1) stop('fs must only have one argument')
  if(length(arg_names) > 2) stop('prop.bkg must have at most two arguments, x and beta')
  if((length(arg_names) == 2) & ((arg_names[1] != 'x') | (arg_names[2] != 'beta'))) stop('Arguments of prop.bkg must be x and beta')
  if((length(arg_names) == 2) & is.null(d_log_bkg)) stop('Must provide the gradient vector d_log_bkg')
  if((length(arg_names) == 2) & is.null(d2_log_bkg)) stop('Must provide the Hessian matrix d2_log_bkg')
  if((length(arg_names) == 2) & is.null(beta_lower)) stop('Must provide the lower limits for beta')
  if((length(arg_names) == 2) & is.null(beta_upper)) stop('Must provide the upper limits for beta')
  if((length(arg_names) == 2) & is.null(beta_start)) stop('Must provide an initial value for beta to obtain MLE')

  if(length(arg_names) == 2){
    arg_names_d_log_bkg <- formalArgs(d_log_bkg)
    if(length(arg_names_d_log_bkg) > 2) stop('d_log_bkg must have at most two arguments')
    if((length(arg_names_d_log_bkg) == 2) & (arg_names_d_log_bkg[2] != 'beta')) stop('The second argument of d_log_bkg must beta')

    if(class(d2_log_bkg)[1] == 'function'){
      arg_names_d2_log_bkg <- formalArgs(d2_log_bkg)
      if(length(arg_names_d2_log_bkg) > 2) stop('d2_log_bkg must have at most two arguments, x and beta')
      if((length(arg_names_d2_log_bkg) == 2) & (arg_names_d2_log_bkg[2] != 'beta')) stop('The second argument of d2_log_bkg must be beta')
    }else if(d2_log_bkg != 0){
      stop("The argument d2_log_bkg must be a function or the numeric value 0")
    }else{
      d2_log_bkg <- matrix(0, p, p)
    }
  }

  # if g needs no parameter estimation:
  if(length(arg_names) == 1){
    prop.bkg <- Vectorize(prop.bkg)
    if(!is.density){
      bkg_mass <- integrate(prop.bkg, l, u)$value
      g <- function(x) prop.bkg(x)/bkg_mass
    }else{
      g <- prop.bkg
    }

    S <- function(x){
      return(fs(x)/g(x) - 1)
    }
    normS <- integrate(function(x){
      (S(x)^2)*g(x)
    }, l,  u)$value |> sqrt()
    S0 <- function(x) {
      S(x)/(normS^2)
    }
    S0_xi <- sapply(xi, S0)

    theta0_hat <- sum(S0_xi*ni)/N; delta0_hat <- sum(S0_xi*mi)/M
    eta_hat <- (theta0_hat - delta0_hat)/(1-delta0_hat)
    test_num <- sqrt(M*N)*(eta_hat - eta0)

    sig_theta0_sq <- sum((S0_xi^2)*ni)/N - theta0_hat^2
    sig_delta0_sq <- sum((S0_xi^2)*mi)/M - delta0_hat^2
    dT_1 <- 1/(1-delta0_hat); dT_2 <- (theta0_hat - 1)/((1-delta0_hat)^2)
    test_denom <- sqrt(M*(dT_1^2)*sig_theta0_sq + N*(dT_2^2)*sig_delta0_sq)
    sig_eta <- test_denom/sqrt(M+N)
    std_err <- sig_eta*sqrt(M+N)/sqrt(M*N)
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
      CI <- c(eta_hat - qnorm(alpha_val)*std_err, Inf)
    }else if(CI_type == 'lower'){
      CI <- c(-Inf, eta_hat + qnorm(alpha_val, lower.tail = FALSE)*std_err)
    }else{
      CI <- eta_hat + c(-1,1)*qnorm(alpha_val/2, lower.tail = FALSE)*std_err
    }
    res <- list(eta_hat = eta_hat, SE = std_err, test.stat = test.stat, p.val = p.val, CI = CI,
                samp.sizes = c(M,N), n_bins = k, alternative = alternative,
                conf.level = conf.level, null.value = eta0, beta.estimate = NULL)
    class(res) <- 'sig_eta_test_bin_wbkg'
    return(res)
  }else if(is.density){
    res <- sig_detect_bin_par_Yden(phys.counts = ni,
                                   bkg.counts = mi,
                                   search_region = search_region,
                                   fs = fs, prop.bkg = prop.bkg,
                                   d_log_bkg = d_log_bkg, d2_log_bkg = d2_log_bkg,
                                   eta0 = eta0,
                                   beta_lower = beta_lower,
                                   beta_upper = beta_upper,
                                   beta_start = beta_start,
                                   conf.level = conf.level,
                                   alternative = alternative,
                                   CI_type = CI_type)
    class(res) <- 'sig_eta_test_bin_wbkg'
    return(res)

  }else{
    res <- sig_detect_bin_par_Nden(phys.counts = ni,
                                   bkg.counts = mi,
                                   search_region = search_region,
                                   fs = fs, prop.bkg = prop.bkg,
                                   d_log_bkg = d_log_bkg, d2_log_bkg = d2_log_bkg,
                                   eta0 = eta0,
                                   beta_lower = beta_lower,
                                   beta_upper = beta_upper,
                                   beta_start = beta_start,
                                   conf.level = conf.level,
                                   alternative = alternative,
                                   CI_type = CI_type)
    class(res) <- 'sig_eta_test_bin_wbkg'
    return(res)
  }
}
