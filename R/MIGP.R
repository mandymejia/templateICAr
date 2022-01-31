#' MIGP
#' 
#' Melodic's Incremental Group PCA
#' 
#' https://doi.org/10.1016/j.neuroimage.2014.07.051
#' 
#' Appendix A
#' 
#' @param dat The files or data to compute PCA for, in list format, with each
#'  subject as a separate list entry.
#' @param datProcFUN Processes each entry of `dat` into a \eqn{T \times V}
#'  matrix (timepoints by locations) which has been column-centered (each 
#'  location's timecourse is mean zero) and variance normalized.
#' @param checkColCentered Check that each subject's data is column centered
#'  after processing it with \code{datProcFUN}? Default: \code{TRUE}.
#' @param nM The size of the "internal" PCA space. Must be larger than the
#'  typical number of timepoints in each individual dataset, \eqn{T}. If
#'  \code{NULL} (default), set to \code{T*2}.
#' @param nP The number of final group PCA components to obtain. If 
#'  \code{NULL} (default), will be set to \code{nM}.
#' @param verbose Occasional updates?
MIGP <- function(dat, datProcFUN, checkColCentered, nM, nP=NULL){
  # Arg checks -----------------------------------------------------------------
  stopifnot(is.list(dat))
  stopifnot(is.function(datProcFUN))

  nN <- length(dat)

  stopifnot(is.numeric(nM))
  stopifnot(length(nM)==1)
  if (nM <= nT) {warning(
    "`nM` should be larger than the ",
    "typical number of timepoints per subject."
  ) }

  if (is.null(nP)) { 
    nP <- nM
  } else {
    stopifnot(is.numeric(nP))
    stopifnot(length(nP)==1)
    stopifnot(nP > 1)
  }

  # First subject --------------------------------------------------------------
  # Read in and process.
  if (verbose) { cat("\tSubject 1: ") }
  dn <- dat[[NN]]
  dn <- datProcFUN(dn)
  nT <- ncol(dn)
  nV <- nrow(dn)
  cat('Number of subjects:            ', nN, "\n")
  cat('Number of data locations:      ', nV, "\n")
  if (verbose) { cat(nT, " timepoints.\n") }

  # Checks
  if (nV > nT) { warning(
    "Data should be TxV after processing, ",
    "but for the first scan there are more rows than columns."
  )}
  if (checkColCentered) {
    if (max(colMeans(dn)) > 1e-8) { 
      stop("First subject's data columns are not demeaned.")
    }
  }

  # Initialize W, the running estimate of the final group-average eigenvectors
  W <- dn # TxV

  # All other subjects ---------------------------------------------------------
  for (nn in seq(2, nN)) {
    if (verbose) { cat("\tSubject ", nn) }
    # Read in and process.
    dn <- dat[[nN]]
    dn <- datProcFUN(dn)
    if (verbose) { cat(ncol(dn), " timepoints.\n") }

    # Checks
    if (nrow(dn) != nV) {
      stop("Subject ", nn, "has ", nrow(dn), " locations (", nV, " expected).")
    }
    if (checkColCentered) {
      if (max(colMeans(dn)) > 1e-8) { 
        stop("Data columns for subject ", nn, " are not demeaned.") 
      }
    }

    # Concatenate
    W <- cbind(W, dn)
    nMm <- min(nM, nrow(W))
    z <- svd(tcrossprod(W), nu=nMm, nv=nMm)
    W <- tcrossprod(diag(z$d), z$v)
  }

  W[seq(nMm),,drop=FALSE]
}