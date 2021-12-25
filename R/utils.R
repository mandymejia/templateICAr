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
#'  matched, or repeated matches occurred, what should happen? Possible values
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
      NULL
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

  INLA_check()
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

#TRY TO REPLACE USE OF THIS WITH REPLACE_XIFTI_DATA
clear_data <- function(x){
  if(!is.null(x$data$cortex_left)) x$data$cortex_left <- matrix(0, nrow(x$data$cortex_left), 1)
  if(!is.null(x$data$cortex_right)) x$data$cortex_right <- matrix(0, nrow(x$data$cortex_right), 1)
  if(!is.null(x$data$subcort)) x$data$subcort <- matrix(0, nrow(x$data$subcort), 1)
  return(x)
}

#' Positive skew?
#' 
#' Does the vector have a positive skew? 
#' 
#' @param x The numeric vector for which to calculate the skew. Can also be a matrix,
#'  in which case the skew of each column will be calculated.
#' @return \code{TRUE} if the skew is positive or zero. \code{FALSE} if the skew is negative.
#' @keywords internal
#' 
#' @importFrom stats median
skew_pos <- function(x){
  x <- as.matrix(x)
  apply(x, 2, median, na.rm=TRUE) <= colMeans(x, na.rm=TRUE)
}

#' Sign match ICA results
#' 
#' Flips all source signal estimates (S) to positive skew
#' 
#' @param x The ICA results with entries \code{S} and \code{M}
#' @return \code{x} but with positive skew source signals
#' @keywords internal
#' 
sign_flip <- function(x){
  stopifnot(is.list(x))
  stopifnot(("S" %in% names(x)) & ("M" %in% names(x)))
  spos <- skew_pos(x$M)
  x$M[,!spos] <- -x$M[,!spos]
  x$S[,!spos] <- x$S[,!spos]
  x
}

#' Center cols
#' 
#' Efficiently center columns of a matrix. (Faster than \code{scale})
#' 
#' @param X The data matrix. Its columns wil lbe centered
#' @return The centered data
#' @keywords internal
colCenter <- function(X) {
  X - rep(colMeans(X), rep.int(nrow(X), ncol(X)))
}

#' Infer fMRI data format
#' 
#' @param BOLD The fMRI data
#' @param verbose Print the format? Default: \code{FALSE}.
#' @return The format: \code{"CIFTI"} file path, \code{"xifti"} object, 
#'  \code{"NIFTI"} file path, \code{"nifti"} object, or \code{"data"}.
#' @keywords internal
infer_BOLD_format <- function(BOLD, verbose=FALSE){

  # Character vector: CIFTI or NIFTI
  if (is.character(BOLD)) {
    format <- ifelse(
      endsWith(BOLD, ".dtseries.nii") | endsWith(BOLD, ".dscalar.nii"),
      "CIFTI", "NIFTI"
    )
    if (length(unique(format)) == 1) {
      format <- format[1]
    } else {
      if (all(endsWith(BOLD, ".nii"))) {
        stop("BOLD format seems to be a mix of CIFTI files and NIFTI files. Use the same format or rename the files.")
      } else {
        stop("BOLD format seems to be a mix of CIFTI files and something else. Use the same format or rename the files.")
      }
    }

  # Non-character vector: xifti, nifti, or data
  } else if (inherits(BOLD[[1]], "xifti")) {
    if (all(lapply(BOLD, inherits, "xifti"))) {
      format <- "xifti"
    } else {
      stop("BOLD format seems to be a mix of `xifti` files and something else. Use the same format for all.")
    }
  } else if (inherits(BOLD[[1]], "nifti")) {
    if (all(lapply(BOLD, inherits, "nifti"))) {
      format <- "nifti"
    } else {
      stop("BOLD format seems to be a mix of `nifti` files and something else. Use the same format for all.")
    }
  } else {
    BOLD_dims <- lapply(BOLD, dim)
    BOLD_dims_lens <- sort(unique(vapply(BOLD_dims, length, 0)))
    if (length(BOLD_dims_lens) > 1) {
      stop("BOLD data have inconsistent dimension lengths. fMRI data should be provided as matrices, not vectors or arrays.")
    } else if (BOLD_dims_lens==4) {
      format <- "nifti" # 4D array: treat as a "nifti"
    } else if (BOLD_dims_lens!=2) {
      stop("BOLD data should be provided as matrices, not vectors or arrays.")
    } else {
      format <- "data"
    }
  }
  if (verbose) { cat("Inferred input format:", format, "\n") }
  format
}

#' Logical argument
#' 
#' Make \code{TRUE} or \code{FALSE}
#' 
#' @param x The argument input
#' @return \code{TRUE} or \code{FALSE}
#' @keywords internal
logical_arg <- function(x, name="x"){
  x <- as.logical(x)
  if (length(x) > 1) { warning("Using the first element of `", name, "`\n"); x <- x[1] }
  x
}

#' Check Q2_max
#' 
#' @param Q2_max,nQ,nT The args
#' @return \code{Q2_max}, clamped to acceptable range of values.
#' @keywords internal
Q2_max_check <- function(Q2_max, nQ, nT){
  if (!is.null(Q2_max)) { 
    if (round(Q2_max) != Q2_max || Q2_max <= 0) {
      stop('`Q2_max` must be `NULL` or a non-negative integer.')
    }
  } else {
    Q2_max <- pmax(round(nT*.50 - nQ), 1)
  }
  
  # This is to avoid the area of the pesel objective function that spikes close 
  #   to rank(X), which often leads to nPC close to rank(X)
  if (Q2_max > round(nT*.75 - nQ)) {
    warning('`Q2_max` too high, setting to 75% of T.')
    Q2_max <- round(nT*.75 - nQ)
  }

  return(Q2_max)
}