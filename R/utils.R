#' Match user inputs to expected values
#'
#' Match each user input to an expected/allowed value. Raise a warning if either
#'  several user inputs match the same expected value, or at least one could not
#'  be matched to any expected value. \code{ciftiTools} uses this function to
#'  match keyword arguments for a function call. Another use is to match
#'  brainstructure labels ("left", "right", or "subcortical").
#'
#' @param user Character vector of user input. These will be matched to
#'  \code{expected} using \code{\link{match.arg}}.
#' @param expected Character vector of expected/allowed values.
#' @param fail_action If any value in \code{user} could not be
#'  matched, or repeated matches occured, what should happen? Possible values
#'  are \code{"stop"} (default; raises an error), \code{"warning"}, and
#'  \code{"nothing"}.
#' @param user_value_label How to refer to the user input in a stop or warning
#'  message. If \code{NULL}, no label is used.
#'
#' @return The matched user inputs.
#'
#' @keywords internal
#'
match_input <- function(
  user, expected,
  fail_action=c("stop", "warning", "message", "nothing"),
  user_value_label=NULL) {

  fail_action <- match.arg(
    fail_action,
    c("stop", "warning", "message", "nothing")
  )
  unrecognized_FUN <- switch(fail_action,
                             stop=stop,
                             warning=warning,
                             message=message,
                             nothing=invisible
  )

  if (!is.null(user_value_label)) {
    user_value_label <- paste0("\"", user_value_label, "\" ")
  }
  msg <- paste0(
    "The user-input values ", user_value_label,
    "did not match their expected values. ",
    "Either several matched the same value, ",
    "or at least one did not match any.\n\n",
    "The user inputs were:\n",
    "\t\"", paste0(user, collapse="\", \""), "\".\n",
    "The expected values were:\n",
    "\t\"", paste0(expected, collapse="\", \""), "\".\n"
  )

  tryCatch(
    {
      matched <- match.arg(user, expected, several.ok=TRUE)
      if (length(matched) != length(user)) { stop() }
      return(matched)
    },
    error = function(e) {
      unrecognized_FUN(msg)
    },
    finally = {
    }
  )

  invisible(NULL)
}



#' Kappa log-likelihood
#' 
#' Compute log-likelihood of kappa given an initial estimate of delta
#'
#' @description Applicable to a single latent field, or multiple latent fields if common smoothness is assumed
#'
#' @param par Vector of length two containing values of log kappa and log residual variance at which to compute log likelihood
#' @param delta Estimate of delta (subject effect or deviation)
#' @param D_diag Diagonal values of D matrix (template standard deviations)
#' @param mesh Object of class "templateICA_mesh" containing the triangular mesh (see \code{\link{make_mesh}})
#' @param C1 For the unit variance case, \eqn{\tau^2 = C1/\kappa^2}, where \eqn{C1 = 1/(4\pi)} when \eqn{\alpha=2}, \eqn{\nu=1}, \eqn{d=2}
#' @param Q Equal to the number of ICs for the common smoothness model, or NULL for the IC-specific smoothness model
#'
#' @return Value of negative log likelihood
#' 
#' @importFrom Matrix bdiag
#' 
#' @keywords internal
#'
loglik_kappa_est <- function(par, delta, D_diag, mesh, C1 = 1/(4*pi), Q=NULL){

  if (!requireNamespace("INLA", quietly = TRUE)) { 
    stop(
      paste0(
        "Package \"INLA\" needed to for spatial modeling.",
        "Please install it at http://www.r-inla.org/download.", 
      ), call. = FALSE
    ) 
  }

  kappa <- exp(par[1]) #log kappa -> kappa
  sigma_sq <- exp(par[2]) #log variance -> variance
  #kappa <- exp(log_kappa)
  #sigma_sq <- exp(log_var)

  Dmat <- Diagonal(length(D_diag), as.vector(D_diag)) #VxV or QVxQV
  delta <- as.vector(delta) #on data locations #length <- V

  #construct indicator matrix of non-data locations in mesh
  Amat <- mesh$A #n_loc x n_mesh
  N <- ncol(mesh$A) #number of mesh locations
  V <- nrow(mesh$A) #number of data locations
  inmesh <- which(colSums(Amat) > 0)
  notinmesh <- setdiff(1:N, inmesh)
  #Imat <- diag(x=1, nrow=N, ncol=N)
  #Amat_c <- Imat[notinmesh,]

  #SPDE matrices
  spde <- mesh$spde
  Fmat <- spde$param.inla$M0
  Gmat <- 1/2*(spde$param.inla$M1 + t(spde$param.inla$M1))
  GFinvG <- spde$param.inla$M2 #this equals G %*% solve(F) %*% G
  Qmat <- C1*(kappa^2 * Fmat + 2 * Gmat + kappa^(-2) * GFinvG)
  Q11 <- Qmat[inmesh,inmesh] # <- Amat %*% Qmat %*% t(Amat)
  if(length(notinmesh) > 0){
    Q12 <- Qmat[inmesh, notinmesh]
    Q21 <- Qmat[notinmesh, inmesh]
    Q22 <- Qmat[notinmesh,notinmesh]
    Q22_inv <- solve(Q22)
    Rinv <- Q11 - (Q12 %*% Q22_inv %*% Q21)
  } else {
    Rinv <- Q11
  }
  cholR <- chol(Rinv) #Rmat <- cholR'cholR, log(det(Rmat)) <- 2*sum(log(diag(cholR)))
  det_Rinv <- 2*sum(log(diag(cholR))) #log determinant
  if(!is.null(Q)) det_Rinv <- Q*det_Rinv

  if(!is.null(Q)) Rinv <- bdiag(rep(list(Rinv), Q))
  W <- Rinv + 1/sigma_sq * (Dmat^2) #W is the matrix K in paper
  cholW <- chol(W) #W <- cholW'cholW
  det_W <- 2*sum(log(diag(cholW))) #log determinant

  #compute determinant part of log-likelihood
  det_sig <- if(is.null(Q)) V*log(sigma_sq) else V*Q*log(sigma_sq)
  det_part <- det_Rinv - det_W - det_sig
  if(abs(det_Rinv) == Inf | abs(det_W) == Inf) {
    stop('negative determinant of precision matrix, returning NA')
    return(NA)
  }

  #compute exponential part of log-likelihood
  D_delta <- Dmat %*% delta
  Winv_D_delta <- INLA::inla.qsolve(Q = W, B=matrix(D_delta, ncol=1), method='solve')
  # mu_post <- 1/sigma_sq * (Dmat %*% Winv_D_delta)
  # Dinv_mupost <- INLA::inla.qsolve(Q = Dmat, B = matrix(mu_post, ncol=1))
  # exp_part1 <- as.numeric(t(Dinv_mupost) %*% Rinv %*% Dinv_mupost)
  # diff <- delta - mu_post
  # exp_part2 <- 1/sigma_sq * sum(diff^2)
  # exp_part <- exp_part1 + exp_part2
  # loglik = det_part - exp_part

  exp_part1 <- as.numeric(1/sigma_sq * sum(delta^2))
  exp_part2 <- as.numeric(1/(sigma_sq^2) * t(D_delta) %*% Winv_D_delta)
  exp_part <- -1*exp_part1 + exp_part2

  loglik <- det_part + exp_part

  return(-1*loglik) #return negative log-likelihood for minimization

}




