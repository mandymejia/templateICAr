#' Create a mask based on vertices that are invalid
#'
#' @param BOLD A \eqn{V \times T} numeric matrix. Each row is a location.
#' @param meanTol,varTol Tolerance for mean and variance of each data location. 
#'  Locations which do not meet these thresholds are masked out of the analysis.
#'  Defaults: \code{-Inf} for \code{meanTol} (ignore), and \code{1e-6} for 
#'  {varTol}.
#' @param verbose Print messages counting how many locations are removed?
#'
#' @importFrom matrixStats rowVars
#' @return A logical vector indicating valid vertices
#'
#' @keywords internal
make_mask <- function(BOLD, meanTol=-Inf, varTol=1e-6, verbose=TRUE){
  stopifnot(is.matrix(BOLD))

  mask_na <- mask_mean <- mask_var <- rep(TRUE, nrow(BOLD))
  # Mark columns with any NA or NaN values for removal.
  mask_na[apply(is.na(BOLD) | is.nan(BOLD), 1, any)] <- FALSE
  # Calculate means and variances of columns, except those with any NA or NaN.
  # Mark columns with mean/var falling under the thresholds for removal.
  mask_mean[mask_na][rowMeans(BOLD[mask_na,,drop=FALSE]) < meanTol] <- FALSE
  if (varTol > 0) {
    mask_var[mask_na][matrixStats::rowVars(BOLD[mask_na,,drop=FALSE]) < varTol] <- FALSE
  }

  # Print counts of locations removed, for each reason.
  if (verbose) {
    warn_part1 <- if (any(!mask_na)) { "additional locations" } else { "locations" }
    if (any(!mask_na)) {
      cat("\t", sum(!mask_na), paste0("locations removed due to NA/NaN values.\n"))
    }
    # Do not include NA locations in count.
    mask_mean2 <- mask_mean | (!mask_na)
    if (any(!mask_mean2)) {
      cat("\t", sum(!mask_mean2), warn_part1, paste0("removed due to low mean.\n"))
    }
    # Do not include NA or low-mean locations in count.
    mask_var2 <- mask_var | (!mask_mean) | (!mask_na)
    if (any(!mask_var2)) {
      cat("\t", sum(!mask_var2), warn_part1, paste0("removed due to low variance.\n"))
    }
  }

  # Return composite mask.
  mask_na & mask_mean & mask_var
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

#' Get FORMAT from format
#' 
#' @param format the file format
#' @return The file FORMAT
#' @keywords internal
#' 
get_FORMAT <- function(format){
  switch(format,
    CIFTI = "CIFTI",
    xifti = "CIFTI",
    GIFTI = "GIFTI",
    gifti = "GIFTI",
    NIFTI = "NIFTI",
    nifti = "NIFTI",
    RDS = "MATRIX",
    data = "MATRIX"
  )
}

#' Check required packages for the data format
#' 
#' @param FORMAT The data FORMAT
#' @return \code{NULL}, invisibly
#' @keywords internal
check_req_ifti_pkg <- function(FORMAT){
  if (FORMAT == "CIFTI") {
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
    }
  }

  if (FORMAT == "GIFTI") {
    if (!requireNamespace("gifti", quietly = TRUE)) {
      stop("Package \"gifti\" needed to work with NIFTI data. Please install it.", call. = FALSE)
    }
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
    }
  }

  if (FORMAT == "NIFTI") {
    if (!requireNamespace("RNifti", quietly = TRUE)) {
      stop("Package \"RNifti\" needed to work with NIFTI data. Please install it.", call. = FALSE)
    }
  }

  invisible(NULL)
}