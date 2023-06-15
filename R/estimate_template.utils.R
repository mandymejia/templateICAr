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

#' Compute theoretical Inverse-Wishart variance of correlation matrix elements
#'
#' @param nu Inverse Wishart degrees of freedom parameter
#' @param p Matrix dimension for IW distribution
#' @param xbar_ij Empirical mean of covariance matrices at element (i,j)
#'
#' @return Theoretical IW variance for correlation element (i,j)
#'
IW_var_cor <- function(nu, p, xbar_ij){
  RHS_numer <- (nu-p+1)*xbar_ij^2 + (nu-p-1)
  RHS_denom <- (nu-p)*(nu-p-3)
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

#' Compute the overall error between empirical and theoretical variance of CORRELATION matrix elements
#'
#' @param nu Inverse Wishart degrees of freedom parameter
#' @param p Matrix dimension for IW distribution
#' @param var Empirical between-subject variances of CORRELATION matrix (upper triangle)
#' @param xbar Empirical mean of CORRELATION matrix (upper triangle)
#' @param M Penalty to assign if theoretical variance is smaller than empirical variance
#'
#' @return Sum of squared difference between the empirical and theoretical IW variance of CORRELATION matrix,
#' but with constraint that theoretical variances must not be smaller than empirical variances
#'
var_sq_err_constrained <- function(nu, p, var, xbar, M=10000){
  K <- length(var)
  var_theo <- rep(NA, K)
  for(k in 1:K){
    var_theo[k] <- IW_var_cor(nu, p, xbar[k])
  }
  penalty <- (var_theo - var)^2
  penalty[var_theo < var] <- M  #constrain so that var_theo >= var
  return(penalty)
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
  for(q1 in 1:(nL-1)){
    for(q2 in (q1+1):nL){
      if(is.na(var_FC[q1,q2])) next()
      #estimate nu by minimizing sq error between the empirical and theoretical IW variance
      nu_opt <- optimize(f=var_sq_err, interval=c(nL+1,nL*100), p=nL, var_ij=var_FC[q1,q2], xbar_ij=mean_FC[q1,q2], xbar_ii=mean_FC[q1,q1], xbar_jj=mean_FC[q2,q2])
      nu_est[q1,q2] <- nu_opt$minimum
    }
  }
  return(nu_est)
}


#' Universally estimate IW dof parameter nu based on method of moments, so
#' that no empirical variance is under-estimated
#'
#' @param var_FC Empirical between-subject variance of covariance matrices (QxQ)
#' @param mean_FC Empirical mean of covariance matrices (QxQ)
#'
#' @importFrom stats optimize
#'
#' @return estimate for nu
#'
estimate_nu <- function(var_FC, mean_FC){

  nL <- nrow(var_FC)
  var_FC_UT <- var_FC[upper.tri(var_FC)]
  mean_FC_UT <- mean_FC[upper.tri(mean_FC)]

  #first, determine the appropriate interval
  interval_LB <- nL+3.1
  interval_UB <- nL*10
  nu_vals <- seq(interval_LB, interval_UB, length.out=1000)
  penalty <- sapply(nu_vals, var_sq_err_constrained, p=nL, var=var_FC_UT, xbar = mean_FC_UT)
  bad_nu <- (colSums(penalty==10000) > 0)

  #if all possible values of nu over-estimate the empirical variance, then return the min val (max var)
  if(sum(!bad_nu) == 0) return(interval_LB)

  #update the UB of the interval to avoid values that under-estimate the empirical var
  interval_UB <- max(nu_vals[!bad_nu]) + (nu_vals[2] - nu_vals[1])
  nu_opt <- optimize(f = function(x){sum(var_sq_err_constrained(x, p=nL, var=var_FC_UT, xbar=mean_FC_UT))},
                     interval=c(interval_LB,interval_UB))
  return(nu_opt$minimum)

}

