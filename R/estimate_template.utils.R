#' Compute theoretical Inverse-Wishart variance of covariance matrix elements
#'
#' @param nu Inverse Wishart degrees of freedom parameter
#' @param p Matrix dimension for IW distribution
#' @param xbar_ij Empirical mean of covariance matrices at element (i,j)
#' @param xbar_ii Empirical mean of covariance matrices at the ith diagonal element
#' @param xbar_jj Empirical mean of covariance matrices at the jth diagonal element
#'
#' @return Theoretical IW variance for covariance element (i,j)
#'
IW_var <- function(nu, p, xbar_ij, xbar_ii, xbar_jj){
  phi_ij <- xbar_ij*(nu-p-1)
  phi_ii <- xbar_ii*(nu-p-1)
  phi_jj <- xbar_jj*(nu-p-1)
  RHS_numer <- (nu-p+1)*phi_ij^2 + (nu-p-1)*phi_ii*phi_jj
  RHS_denom <- (nu-p)*((nu-p-1)^2)*(nu-p-3)
  return(RHS_numer/RHS_denom)
}

#' Compute the error between empirical and theoretical variance of covariance matrix elements
#'
#' @param nu Inverse Wishart degrees of freedom parameter
#' @param p Matrix dimension for IW distribution
#' @param var_ij Empirical between-subject variance of covariance matrices at element (i,j)
#' @param xbar_ij Empirical mean of covariance matrices at element (i,j)
#' @param xbar_ii Empirical mean of covariance matrices at the ith diagonal element
#' @param xbar_jj Empirical mean of covariance matrices at the jth diagonal element
#'
#' @return Squared difference between the empirical and theoretical IW variance of covariance matrices at element (i,j)
#'
var_sq_err <- function(nu, p, var_ij, xbar_ij, xbar_ii, xbar_jj){
  var_ij_theo <- IW_var(nu, p, xbar_ij, xbar_ii, xbar_jj)
  sq_diff <- (var_ij - var_ij_theo)^2
  return(sq_diff)
}

#' Estimate IW dof parameter nu based on method of moments
#'
#' @param var_FC Empirical between-subject variance of covariance matrices (QxQ)
#' @param mean_FC Empirical mean of covariance matrices (QxQ)
#'
#' @importFrom stats optimize
#'
#' @return QxQ matrix of estimates for nu
#'
estimate_nu_matrix <- function(var_FC, mean_FC){
  nL <- nrow(var_FC)
  nu_est <- matrix(NA, nL, nL)
  for(q1 in 1:nL){
    for(q2 in q1:nL){
      if(is.na(var_FC[q1,q2])) next()
      #estimate nu by minimizing sq error between the empirical and theoretical IW variance
      nu_opt <- optimize(f=var_sq_err, interval=c(nL+1,nL*10), p=nL, var_ij=var_FC[q1,q2], xbar_ij=mean_FC[q1,q2], xbar_ii=mean_FC[q1,q1], xbar_jj=mean_FC[q2,q2])
      nu_est[q1,q2] <- nu_opt$minimum
    }
  }
  return(nu_est)
}
