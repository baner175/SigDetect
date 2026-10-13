#' Calculates the signal region with mass `1-eps`
#' @param Fs the function for the cumulative distribution function of the signal.
#' @param mu_sig location of the signal distribution
#' @param search_region a numeric vector of two numbers with the first and the second number being the lower and upper limits of the search region respectively.
#' @param eps area outside the signal region to consider.
#' @return a numeric vector `c(M_lower, M_upper)` the lower and upper limits of the signal region with mass `1-eps`.
#'
#' @export
signal_region <- function(Fs, mu_sig,
                          search_region, eps = 0.001){
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
  return(c(M_lower, M_upper))
}
