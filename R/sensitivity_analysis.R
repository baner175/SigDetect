
sensitivity_analysis_unbinned <- function(bench.bkg,
                                          search_region,
                                          Fs, mu_sig,
                                          is.density = TRUE,
                                          phys.samp = NULL,
                                          eps = 0.001,
                                          lambda,
                                          mu1 = NULL, mu2 = NULL, sig0,
                                          alpha_start = NULL,
                                          alpha_lower = NULL,
                                          alpha_upper = NULL,
                                          scale_by,
                                          ...){
  q_names <- formalArgs(bench.bkg)
  if(length(lambda)>5) stop("At most five values of lambda can be used")
  if(!is.null(phys.samp)){
    if(is.density){
      neg_loglik <- function(a){
        qi <- sapply(phys.samp, bench.bkg, alpha = a)
        return(-sum(log(qi)))
      }
      alpha_hat <- nlminb(start = alpha_start,
                          objective = neg_loglik,
                          lower = alpha_lower,
                          upper = alpha_upper)$par
      q <- function(x) bench.bkg(x, alpha = alpha_hat)
    }else{
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
    }
  }else if(length(q_names) > 1){
    stop('When Physics data is not available, bench.bkg must only have one argument')
  }else{
    if(is.density){
      q <- bench.bkg
    }else{
      q_mass <- integrate(bench.bkg, l, u)$value
      q <- function(x) bench.bkg(x)/q_mass
    }
  }


  l <- search_region[1]; u <- search_region[2]
  find_d <- function(d)
  {
    pl <- Fs(mu_sig-d)
    pu <- Fs(mu_sig+d)
    return(pu-pl-1+eps)
  }
  sol <- uniroot(find_d, lower = 0, upper = min(mu_sig - l,u - mu_sig))
  r <- sol$root
  M_lower <- mu_sig - r
  M_upper <- mu_sig + r

  if(is.null(mu1)) mu1 <- (mu_sig+M_lower)/2
  if(is.null(mu2)) mu2 <- (mu_sig+M_upper)/2

  if(!is.null(phys.samp)){
    hist(phys.samp, breaks = 30, probability = TRUE,
         main = 'Histogram of the physics sample',
         xlab = 'x', ylab = 'Density')
    curve(q, lwd = 4, lty = 1, col = 'black', add = TRUE)
  }else{
    curve(q(x)*scale_by, lwd = 4, lty = 1, col = 'black', from = l, to = u,
          ylab = 'Density', xlab = 'x')
  }

  mycols <- c('blue', 'brown', 'darkgreen', 'orange', 'red')
  my_lty <- c(3,2,6,5,4)

  for(j in 1:length(lambda))
  {
    curve(((1-2*lambda[j])*q(x)+
             lambda[j]*(truncdist::dtrunc(x, spec = 'norm', a = l, b = u,
                                          mean = mu1, sd = sig0) +
                          truncdist::dtrunc(x, spec = 'norm', a = l, b = u,
                                            mean = mu2, sd = sig0)))*scale_by,
          add = TRUE, lwd = 4,
          col = mycols[j],
          lty = my_lty[j])
  }

  abline(v = c(M_lower, M_upper), col = 'grey', lty = 2, lwd = 4)
  polygon(
    x = c(M_lower, M_upper, M_upper, M_lower),
    y = c(0, 0, 1e10, 1e10),
    col = rgb(0, 1, 0, 0.15),
    border = NA)

  legend('topright',
         col = c('black', mycols),
         lty = c(1,my_lty), bty = 'n', lwd = 4,
         legend = c(
           expression(lambda == 0),
           sapply(lambda, function(x) bquote(lambda == .(x)))),
         cex = 2,
         x.intersp = 0.5,
         seg.len = 2,
         y.intersp = 1)

}


sensitivity_analysis_binned <- function(bench.bkg,
                                        search_region,
                                        Fs, mu_sig,
                                        is.density = TRUE,
                                        phys.counts = NULL,
                                        eps = 0.001,
                                        lambda,
                                        mu1 = NULL, mu2 = NULL, sig0,
                                        alpha_start = NULL,
                                        alpha_lower = NULL,
                                        alpha_upper = NULL,
                                        scale_by,
                                        ...){
  q_names <- formalArgs(bench.bkg)
  l <- search_region[1]; u <- search_region[2]
  if(length(lambda)>5) stop("At most five values of lambda can be used")
  if(!is.null(phys.counts)){
    ni <- phys.counts
    k <- length(phys.counts)
    bin_ends <- seq(l, u, length.out = k+1)
    ni <- phys.counts
    if(is.density){
      neg_loglik <- function(a){
        qi <- sapply(1:k, function(i){
          integrate(bench.bkg, bin_ends[i], bin_ends[i+1], alpha = a)$value
        })
        return(-sum(ni*log(qi)))
      }
      alpha_hat <- nlminb(start = alpha_start,
                          objective = neg_loglik,
                          lower = alpha_lower,
                          upper = alpha_upper)$par
      q <- function(x) bench.bkg(x, alpha = alpha_hat)
    }else{
      neg_loglik <- function(a){
        h_mass <- integrate(bench.bkg, l, u, alpha = a)$value
        qi <- sapply(1:k, function(i){
          integrate(bench.bkg, bin_ends[i], bin_ends[i+1], alpha = a)$value/h_mass
        })
        return(-sum(ni*log(qi)))
      }
      alpha_hat <- nlminb(start = alpha_start,
                          objective = neg_loglik,
                          lower = alpha_lower,
                          upper = alpha_upper)$par
      h_mass <- integrate(bench.bkg, l, u, alpha = alpha_hat)$value
      h <- function(x) bench.bkg(x, alpha = alpha_hat)
      q <- function(x) bench.bkg(x, alpha = alpha_hat)/h_mass
    }
  }else if(length(q_names) > 1){
    stop('When Physics data is not available, bench.bkg must only have one argument')
  }else{
    if(is.density){
      q <- bench.bkg
    }else{
      q_mass <- integrate(bench.bkg, l, u)$value
      q <- function(x) bench.bkg(x)/q_mass
    }
  }

  find_d <- function(d)
  {
    pl <- Fs(mu_sig-d)
    pu <- Fs(mu_sig+d)
    return(pu-pl-1+eps)
  }
  sol <- uniroot(find_d, lower = 0, upper = min(mu_sig - l,u - mu_sig))
  r <- sol$root
  M_lower <- mu_sig - r
  M_upper <- mu_sig + r

  if(is.null(mu1)) mu1 <- (mu_sig+M_lower)/2
  if(is.null(mu2)) mu2 <- (mu_sig+M_upper)/2
  mycols <- c('blue', 'brown', 'darkgreen', 'orange', 'red')
  my_lty <- c(3,2,6,5,4)
  if(!is.null(phys.counts)){
    N <- sum(ni)
    xi <- (bin_ends[1:k] + bin_ends[2:(k+1)])/2
    plot(x = xi,  y = ni,
         main = '',
         xlab = 'x', ylab = 'Counts')

    curve(N*q(x)*(u-l)/(k+1), lwd = 4, lty = 1, col = 'black', add = TRUE)
    # rect(
    #   xleft = M_lower,
    #   xright = M_upper,
    #   ybottom = -1e10,
    #   ytop = 1e10,
    #   col = rgb(0, 1, 0, 0.15),
    #   border = NA
    # )
    for(j in 1:length(lambda))
    {
      curve(((1-2*lambda[j])*q(x)+
               lambda[j]*(truncdist::dtrunc(x, spec = 'norm', a = l, b = u,
                                            mean = mu1, sd = sig0) +
                            truncdist::dtrunc(x, spec = 'norm', a = l, b = u,
                                              mean = mu2, sd = sig0)))*N*(u-l)/(k+1),
            add = TRUE, lwd = 4,
            col = mycols[j],
            lty = my_lty[j])
    }

  }else{
    curve(q(x)*scale_by, lwd = 4, lty = 1, col = 'black', from = l, to = u,
          ylab = 'Density', xlab = 'x')
    # rect(
    #   xleft = M_lower,
    #   xright = M_upper,
    #   ybottom = -1e10,
    #   ytop = 1e10,
    #   col = rgb(0, 1, 0, 0.15),
    #   border = NA
    # )
    for(j in 1:length(lambda))
    {
      curve(((1-2*lambda[j])*q(x)+
               lambda[j]*(truncdist::dtrunc(x, spec = 'norm', a = l, b = u,
                                            mean = mu1, sd = sig0) +
                            truncdist::dtrunc(x, spec = 'norm', a = l, b = u,
                                              mean = mu2, sd = sig0)))*scale_by,
            add = TRUE, lwd = 4,
            col = mycols[j],
            lty = my_lty[j])
    }
  }

  abline(v = c(M_lower, M_upper), col = 'grey', lty = 2, lwd = 4)
  polygon(
    x = c(M_lower, M_upper, M_upper, M_lower),
    y = c(0, 0, 1e10, 1e10),
    col = rgb(0, 1, 0, 0.15),
    border = NA)
  legend('topright',
         col = c('black', mycols),
         lty = c(1,my_lty), bty = 'n', lwd = 4,
         legend = c(
           expression(lambda == 0),
           sapply(lambda, function(x) bquote(lambda == .(x)))),
         cex = 2,
         x.intersp = 0.5,
         seg.len = 2,
         y.intersp = 1)

}
#' Performs sensitivity analysis using the benchmark background model
#' @param bench.bkga function specifying the benchmark background model with arguments `x` and `alpha`. The function is not required to be a density.
#' @param search_region a numeric vector of two numbers with the first and the second number being the lower and upper limits of the search region respectively.
#' @param Fs the function for the cumulative distribution function of the signal.
#' @param mu_sig location of the signal distribution
#' @param is.density a logical indicating whether `bench.bkg` is a density.
#' @param phys.data a numeric vector containing the physics data (binned or unbinned)
#' @param is.binned a logical indicating whether `phys.data` is binned.
#' @param eps area outside the signal region to consider.
#' @param lambda a numeric vector containing the the sizes of the diffused bump added to `bench.bkg`. The values must be between 0 and 0.5. The vector should not contain more that five values.
#' @param mu1 a numeric value for the location of the first Gaussian component in the diffused bump.Default is `NULL`. When `NULL`, `mu1` is the average of the lower limit of the signal region and `mu_sig`.
#' @param mu2 a numeric value for the location of the second Gaussian component in the diffused bump.Default is `NULL`. When `NULL`, `mu2` is the average of the lower limit of the signal region and `mu_sig`.
#' @param sig0 a numeric value for the common scale parameter of the Gaussian components in the diffused bump.
#' @param alpha_lower a numeric vector or a numeric value for the lower limit of `alpha`, the parameter in `bench.bkg`.
#' @param alpha_upper a numeric vector or a numeric value for the upper limit of `alpha`, the parameter in `bench.bkg`.
#' @param alpha_start a numeric vector or a numeric value for the initial value of `alpha` to obtains its MLE.
#' @param scale_by a numeric value to scale the densities in the generated plot.
#'
#' @export

sensitivity_analysis <- function(bench.bkg,
                                 search_region,
                                 Fs, mu_sig,
                                 is.density = TRUE,
                                 phys.data = NULL,
                                 is.binned = FALSE,
                                 eps = 0.001,
                                 lambda,
                                 mu1 = NULL, mu2 = NULL, sig0,
                                 alpha_start = NULL,
                                 alpha_lower = NULL,
                                 alpha_upper = NULL,
                                 scale_by = 1,
                                 ...){
  n_lamb_out <- sum((lambda < 0) | (lambda > 0.5))
  if(n_lamb_out > 0) stop('lambda values must be between 0 and 0.5')
  if(is.binned){
    sensitivity_analysis_binned(bench.bkg = bench.bkg,
                                search_region = search_region,
                                Fs = Fs, mu_sig = mu_sig,
                                is.density = is.density,
                                phys.counts = phys.data,
                                eps = eps,
                                lambda = lambda,
                                mu1 = mu1, mu2 = mu2, sig0 = sig0,
                                alpha_start = alpha_start,
                                alpha_lower = alpha_lower,
                                alpha_upper = alpha_upper,
                                scale_by = scale_by,
                                ...)
  }else{
    sensitivity_analysis_unbinned(bench.bkg = bench.bkg,
                                  search_region = search_region,
                                  Fs = Fs, mu_sig = mu_sig,
                                  is.density = is.density,
                                  phys.samp = phys.data,
                                  eps = eps,
                                  lambda,
                                  lambda = lambda,
                                  mu1 = mu1, mu2 = mu2, sig0 = sig0,
                                  alpha_start = alpha_start,
                                  alpha_lower = alpha_lower,
                                  alpha_upper = alpha_upper,
                                  scale_by = scale_by,,
                                  ...)
  }
}
