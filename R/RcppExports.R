# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Update parameters (M-step of the EM)
#'
#' @param template_mean a matrix with dimensions V x Q giving the mean value
#'   of the independent components
#' @param template_var a matrix with dimensions V x Q giving the variance
#'   of the independent components
#' @param template_FC a list with two elements: psi, a Q x Q matrix, and nu, a
#'   scalar. These two values are the parameters of the Wishart prior on G
#' @param G a Q x Q matrix of the prior covariance of A
#' @param prior_params a length 2 vector with the prior parameters for tau_v
#' @param BOLD a V x T matrix of BOLD values
#' @param Y_sq_sum a length V vector with the sum of squared BOLD values at
#'   each data location
#' @param post_sums a list of posterior summary statistics including the named
#'   summaries \code{AS_sq_sum}, \code{yAS_sum}, \code{A_sum}, and \code{AtA_sum}
#' @param sigma2_alpha a scalar multiplier for the prior variance of alpha
#' @param verbose a boolean. Should messages be generated and output?
#' @return A list with quantities tau_sq, alpha, and G
#' @export
UpdateTheta_FCtemplateICAcpp <- function(template_mean, template_var, template_FC, G, prior_params, BOLD, Y_sq_sum, post_sums, sigma2_alpha, verbose) {
    .Call(`_templateICAr_UpdateTheta_FCtemplateICAcpp`, template_mean, template_var, template_FC, G, prior_params, BOLD, Y_sq_sum, post_sums, sigma2_alpha, verbose)
}

#' Use a Gibbs sampler for the A and S variables (E-step of the EM)
#'
#' @param nsamp the number of posterior samples to output after burn-in
#' @param nburn the number of posterior samples to throw away before saving
#' @param template_mean a matrix with dimensions V x Q giving the mean value
#'   of the independent components
#' @param template_var a matrix with dimensions V x Q giving the variance
#'   of the independent components
#' @param S a matrix with dimensions V x Q of subject independent components
#' @param G a Q x Q matrix of the prior covariance of A
#' @param tau_v a length V vector with noise variance for each data location
#' @param Y a matrix with dimensions V x T of observed BOLD data
#' @param alpha a length Q vector of the prior mean of all rows of A
#' @param final a boolean. Should posterior samples be returned instead of
#'   summary measures?
#' @param return_samp a boolean. Should posterior samples be returned?
#' @return List with estimates for A, S, and possibly other quantities
#' @export
Gibbs_AS_posteriorCPP <- function(nsamp, nburn, template_mean, template_var, S, G, tau_v, Y, alpha, final, return_samp) {
    .Call(`_templateICAr_Gibbs_AS_posteriorCPP`, nsamp, nburn, template_mean, template_var, S, G, tau_v, Y, alpha, final, return_samp)
}

