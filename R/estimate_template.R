#' Estimate template from DR estimates
#'
#' @param DR1,DR2 the test and retest lists of dual regression estimates
#' @param var_method \code{"unbiased"} (default) or \code{"non-negative"}
#'
#' @return List of two elements: the mean and variance templates
#' @keywords internal
#' @importFrom abind abind
estimate_template_from_DR <- function(
  DR1, DR2, var_method=c("unbiased", "non-negative")){

  # Check arguments.
  stopifnot(length(dim(DR1)) == length(dim(DR2)))
  stopifnot(all(dim(DR1) == dim(DR2)))
  N <- dim(DR1)[1]
  var_method <- match.arg(var_method, c("unbiased", "non-negative"))

  template <- list(mean=NULL, var=NULL)

  # Mean.
  template$mean <- t(colMeans(DR1 + DR2, na.rm=TRUE) / 2)

  # Variance.
  # Unbiased: hat(sigmasq_btwn)
  # Non-negative: MSB_div2 === hat(sigmasq_btwn) + hat(sigmasq_noise) / k
  SSB <- 2 * colSums(((DR1 + DR2)/2 - rep(t(template$mean), each=N))^2, na.rm=TRUE)
  MSB_div2 <- t(SSB / (N-1)) / 2
  if (var_method == "unbiased") {
    # Fastest method.
    var_noise <- t( (1/2) * apply(DR1 - DR2, c(2,3), var, na.rm=TRUE) )
    template$var <- MSB_div2 - var_noise/2

    # # Previous, equivalent calculation.
    # var_tot1 <- apply(DR1, c(2,3), var, na.rm=TRUE)
    # var_tot2 <- apply(DR2, c(2,3), var, na.rm=TRUE)
    # var_tot <- t((var_tot1 + var_tot2)/2)
    # # noise (within-subject) variance
    # DR_diff <- DR1 - DR2;
    # var_noise <- t((1/2)*apply(DR_diff, c(2,3), var, na.rm=TRUE))
    # # signal (between-subject) variance
    # template$var <- var_tot - var_noise

    # # Another equivalent calculation.
    # template$var <- t(apply(
    #   abind::abind(DR1, DR2, along=1),
    #   seq(2, 3),
    #   function(q){ cov(q[seq(N)], q[seq(N+1, 2*N)], use="complete.obs") }
    # ))

    # Make negative estimates equal to zero.
    template$var[template$var < 0] <- 0

  } else {
    template$var <- MSB_div2
  }

  template
}
